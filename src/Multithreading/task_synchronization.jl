"""
    spawn_tasks_safe_wait(f, num_tasks::Int)

Run `f(taskid::Int, abort::Threads.Atomic{Bool})` as `num_tasks` spawned tasks, and wait for
all of them to actually finish before returning or throwing - whether that is because they
succeeded, failed, or observed `abort[]` and stopped early. `f` must check `abort[]` itself;
this function does not interrupt it.

On failure, rethrows the first task's exception (by actual failure order, see below) as a
plain `TaskFailedException`, never a `CompositeException` (unlike `@sync`).

**Developer notes:** the CAS recording a failure happens *before* `abort[]` is set, so a
sibling that observes `abort[]` and only then fails itself can never win that CAS against the
failure it just observed. A task's own failure takes precedence over an asynchronous
exception observed while waiting (see [`drain_tasks`](@ref)), even if the latter happened
first - the semantics are "this call failed", with asynchronous exceptions only reported when
nothing failed on its own.
"""
function spawn_tasks_safe_wait(f, num_tasks::Int)
    abort = Threads.Atomic{Bool}(false)
    first_failed = Threads.Atomic{Int}(0)
    tasks = Vector{Task}(undef, num_tasks)
    for taskid in 1:num_tasks
        tasks[taskid] = Threads.@spawn begin
            try
                f(taskid, abort)
            catch
                Threads.atomic_cas!(first_failed, 0, taskid) # compare-and-set (if not set)
                abort[] = true
                rethrow()
            end
        end
    end
    async_exceptions = drain_tasks(tasks)
    n_async = length(async_exceptions)

    failed_taskid = first_failed[]
    if failed_taskid != 0
        n_async == 0 || @info """
            $n_async asynchronous exception(s) (e.g. a Ctrl-C interrupt) arrived while
            waiting for tasks to finish; a task also failed, which takes precedence
            """
        throw(TaskFailedException(tasks[failed_taskid]))
    elseif !isempty(async_exceptions)
        if n_async > 1
            types = join(string.(unique(typeof.(async_exceptions))), ", ")
            @info """
                $n_async asynchronous exceptions arrived while waiting for tasks to finish.
                The different types of exceptions where {$types},
                and first exception is propagated (thrown).
                """
        end
        throw(first(async_exceptions))
    end
    return nothing
end

"""
    drain_tasks(tasks::AbstractVector{<:Task}) -> Vector{Exception}

Wait for every task in `tasks` to actually finish, returning any exceptions that were
delivered *asynchronously to this calling task* while doing so (e.g. a `Ctrl-C` interrupt),
in the order they were observed.

**Developer notes:** `wait(t::Task)` throws `TaskFailedException(t)` if `t` failed.
Any other exception caught from `wait(t)` is asynchronous (did not come from `t`),
and `t` may still be running - this includes a `TaskFailedException` wrapping some
*other* task: that can be delivered asynchronously too (e.g. forwarded by external
supervising code), so checking the exception's own `.task` field, not just its type,
is what actually proves `t` has failed. `Base.waitall` will be killed by an
asynchronous exception and may leave tasks running.
"""
function drain_tasks(tasks::AbstractVector{<:Task})
    async_exceptions = Exception[]
    for t in tasks
        while true
            try
                wait(t; throw=false)
            catch e
                push!(async_exceptions, e)
                continue
            end
            break
        end
    end
    return async_exceptions
end
