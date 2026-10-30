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

An `InterruptException` delivered asynchronously to the *calling* task (e.g. Ctrl-C) while it
is waiting on one of the spawned tasks is distinguished from that task failing on its own:
`wait` can throw for either reason, and only an `InterruptException` thrown directly into the
waiting call (as opposed to a `TaskFailedException` wrapping some other error from the task
itself) is retried, on the same task, until it actually finishes. The interrupt is remembered
and only re-thrown (to the caller of `sync_with_abort`) once every task has been drained this
way - so it still reaches the caller instead of being silently swallowed, but never at the
cost of returning while a task is still running.
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
                # Order matters: the CAS must happen *before* `abort[] = true`, not after.
                # Any task observing `abort[]` can then be certain that whichever task set it
                # had already completed its own CAS attempt first (that's a same-task program-
                # order guarantee, unlike the relative timing of two different tasks) - so a
                # sibling that waits on `abort[]` before failing itself can never win the CAS
                # against a failure that was already underway when it started waiting.
                Threads.atomic_cas!(first_failed, 0, taskid)
                abort[] = true
                rethrow()
            end
        end
    end
    interrupted = false
    foreach(tasks) do t
        while true
            try
                wait(t)
                break
            catch e
                if e isa InterruptException
                    interrupted = true
                    continue
                end
                break
            end
        end
    end
    failed_taskid = first_failed[]
    if failed_taskid != 0
        throw(TaskFailedException(tasks[failed_taskid]))
    elseif interrupted
        throw(InterruptException())
    end
    return nothing
end
