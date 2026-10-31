"""
    sync_with_abort(f, num_tasks::Int)

Run `f(taskid::Int, abort::Threads.Atomic{Bool})` concurrently as `num_tasks` spawned tasks
(`taskid in 1:num_tasks`), sharing the single `abort` flag between them.

This exists as a restricted, purpose-built replacement for `Base.Experimental.@sync`
(fixing the specific defect that makes it unsafe for a use case like this one, where the
tasks share long-lived, task-local state across repeated calls - see
[`TaskLocals`](@ref)/`ThreadedDomainBuffer`): `sync_with_abort` always waits for *every* task
to finish - successfully, having stopped early because `abort[]` became `true`, or having
itself thrown - before this function returns or propagates an exception. `f` is responsible
for checking `abort[]` itself and returning early when it is set; `sync_with_abort` does not
interrupt `f`.

If one or more invocations of `f` throw, every task is still waited for (`abort[]` is set so
well-behaved siblings stop early instead of continuing wasted work), and then one exception is
rethrown, wrapped in a plain `TaskFailedException` (never a `CompositeException`, unlike plain
`@sync`): the first one *recorded* by a `compare_and_set!` race between the tasks' own `catch`
blocks (each does the CAS *before* setting `abort[]`, so a sibling that observes `abort[]` and
only then fails itself can never win this race against the failure it just observed - it is
only a race between failures with no such causal relationship between them). In practice this
is the first task to actually fail, but note that if two tasks fail near-simultaneously with
neither having observed the other's `abort[]` write yet, this is decided by which `catch`
block's CAS runs first, not a hardware-level global clock ordering. A task's own failure always takes
precedence over a remembered caller interrupt (see below), even if the interrupt happened
"first": the semantics are "this call failed", with the interrupt as a secondary detail about
how the wait went.

Any exception delivered *asynchronously to the calling task* (e.g. an `InterruptException`
from Ctrl-C, or some other exception thrown in by a third party) while it is waiting on one
of the spawned tasks is distinguished from that task failing on its own: `wait(t::Task)`
throws `TaskFailedException` for exactly one reason - `t` itself failed - so any other
exception caught from it cannot say anything about `t` (which may still be running,
completely unaffected) and is retried, on the same task, until it actually finishes (see
[`drain_tasks`](@ref)). Such exceptions are remembered and only re-thrown (to the caller of
`sync_with_abort`) once every task has been drained this way - so they still reach the caller
instead of being silently swallowed, but never at the cost of returning while a task is still
running.
"""
function sync_with_abort(f, num_tasks::Int)
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

**Developer notes:** `wait(t::Task)` throws `TaskFailedException` if `t` failed. 
Any other exception caught from `wait(t)` is asynchronous (did not come from `t`), 
and `t` may still be running. `waitall` will be killed by an asynchronous exception
and may leave tasks running.
"""
function drain_tasks(tasks::AbstractVector{<:Task})
    async_exceptions = Exception[]
    for t in tasks
        while true
            try
                wait(t)
            catch e
                if !isa(e, TaskFailedException)
                    push!(async_exceptions, e)
                    continue
                end
            end
            break
        end
    end
    return async_exceptions
end
