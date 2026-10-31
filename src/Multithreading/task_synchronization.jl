"""
    spawn_tasks_safe_wait(f, num_tasks::Int)

Run `f(taskid::Int, abort::Threads.Atomic{Bool})` as `num_tasks` spawned tasks, and wait for
all of them to actually finish before returning or throwing - whether that is because they
succeeded, failed, or observed `abort[]` and stopped early. `f` must check `abort[]` itself;
this function does not interrupt it.

### Difference to `@sync` + `@spawn`
`spawn_tasks_safe_wait(f, num_tasks)` is very similar to the following code,
```julia
abort = Threads.Atomic{Bool}(false)
@sync begin
    for taskid in 1:num_tasks
        Threads.@spawn begin
            f(taskid, abort)
        end
    end
end
```
But there are a few key differences
* If killed by an asynchronous exception, the function will still wait for all tasks to finish.
  `@sync` will only do this for exceptions thrown from one of the spawned tasks, and
  `Base.Experimental.@sync` will never wait.
* The returned error if a task throws is `TaskFailedException`, as opposed to the `CompositeException`
  thrown with `@sync`. (`Base.Experimental.@sync` throws a `TaskFailedException`)
* `abort[]` becomes `true` when any task fails, allowing the other tasks to end early.
  This holds even if the caller task throws (e.g. if Ctrl-C arrives while waiting for tasks to finish)

### Developer notes
`compare-and-set` of `first_failed` happens *before* `abort[] = true` is set, so a
sibling failing after seeing `abort[] == true` cannot set `first_failed`.

A task's own failure takes precedence over an asynchronous exception observed while waiting
(see [`drain_tasks`](@ref)), even if the latter happened first.

`disable_sigint` ensures that a Ctrl-C event cannot land *between* spawning a task and recording 
it in `tasks` (losing its reference so `drain_tasks` doesn't see it).
Exception raised while spawning is treated as asynchronous exception: 
`abort[]=true` and all spawned tasks are drained before it is reported.
"""
function spawn_tasks_safe_wait(f, num_tasks::Int)
    abort = Threads.Atomic{Bool}(false)
    first_failed = Threads.Atomic{Int}(0)
    tasks = Task[]
    async_exceptions = Exception[]
    try
        disable_sigint() do
            for taskid in 1:num_tasks
                push!(tasks, Threads.@spawn begin
                    try
                        f(taskid, abort)
                    catch
                        Threads.atomic_cas!(first_failed, 0, taskid) # compare-and-set (if not set)
                        abort[] = true
                        rethrow()
                    end
                end)
            end
        end
    catch e
        abort[] = true
        push!(async_exceptions, e)
    end
    append!(async_exceptions, drain_tasks(tasks, abort))
    n_async = length(async_exceptions)

    failed_taskid = first_failed[]
    if failed_taskid != 0
        n_async == 0 || @info """
            $n_async asynchronous exception(s) (e.g. a Ctrl-C interrupt) arrived while
            spawning or waiting for tasks to finish; a task also failed, which takes
            precedence
            """
        throw(TaskFailedException(tasks[failed_taskid]))
    elseif !isempty(async_exceptions)
        if n_async > 1
            types = join(string.(unique(typeof.(async_exceptions))), ", ")
            @info """
                $n_async asynchronous exceptions arrived while spawning or waiting for tasks
                to finish. The different types of exceptions where {$types},
                and first exception is propagated (thrown).
                """
        end
        throw(first(async_exceptions))
    end
    return nothing
end

"""
    drain_tasks(tasks::AbstractVector{<:Task}, abort::Threads.Atomic{Bool}) -> Vector{Exception}

Wait for every task in `tasks` to actually finish, returning any exceptions that were
delivered *asynchronously to this calling task* while doing so (e.g. a `Ctrl-C` interrupt),
in the order they were observed. If any such exceptions are detected, sets `abort[] = true`,
allowing tasks that observe this to finish early.

`current_task() ∈ tasks` is not allowed, and will throw a `ConcurrencyViolationError`
after draining the other tasks in `tasks`.

### Developer notes
`wait(t; throw = false)` blocks until `t` is actually done -
succeeded or failed, it doesn't matter which, that isn't this function's concern - and,
per its own contract, never throws because of `t`'s *own* outcome. So *any* exception
caught here is necessarily asynchronous (did not come from `t`, which may still be
running), with no need to inspect the exception's type or identity to tell the two cases
apart. `Base.waitall` will be killed by an asynchronous exception and may leave tasks
running.
"""
function drain_tasks(tasks::AbstractVector{<:Task}, abort::Threads.Atomic{Bool})
    async_exceptions = Exception[]
    for t in tasks
        t === current_task() && continue # Just skip it, cannot wait for itself
        while true
            try
                wait(t; throw=false)
            catch e
                abort[] = true
                push!(async_exceptions, e)
                continue
            end
            break
        end
    end
    if any(t === current_task() for t in tasks)
        throw(ConcurrencyViolationError(
            """`drain_tasks` cannot be called with the current task in `tasks`:
            The calling task cannot be drained. All other tasks have been drained"""))
    end
    return async_exceptions
end
