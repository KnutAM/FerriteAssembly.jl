# Fast, isolated tests for `FerriteAssembly.spawn_tasks_safe_wait` (src/Multithreading/
# task_synchronization.jl), the primitive extracted from `work_domain_threaded!` that
# replaces `Base.Experimental.@sync` for this use case. No grid/CellValues/DomainBuffer
# machinery is needed here: these tests exercise the primitive directly with trivial
# closures, so they're quick and - unlike driving the same guarantees through a full
# threaded assembly - not sensitive to scheduling noise under load.
#
# Ordering between tasks is established with explicit handshakes (a task waits on a flag
# another task sets), never by guessing relative `sleep` durations: a wrong guess could make
# a test pass vacuously (e.g. a "slow" sibling that happens to finish before the "fast"
# failure) instead of actually exercising the intended race.
@testset "spawn_tasks_safe_wait" begin
    function wait_for(cond; timeout = 10.0)
        t0 = time()
        while !cond()
            time() - t0 > timeout && return false
            sleep(0.001) # see test/threading_errors.jl for why not `yield()`
        end
        return true
    end

    @testset "all tasks succeed" begin
        calls = Threads.Atomic{Int}(0)
        seen = fill(false, 4)
        lk = ReentrantLock()
        FerriteAssembly.spawn_tasks_safe_wait(4) do taskid, abort
            Threads.atomic_add!(calls, 1)
            lock(lk) do
                seen[taskid] = true
            end
        end
        @test calls[] == 4
        @test all(seen)
    end

    @testset "one failure: siblings observe abort, all are waited for" begin
        done = Threads.Atomic{Int}(0)
        abort_observed = fill(false, 3)
        e = nothing
        try
            FerriteAssembly.spawn_tasks_safe_wait(4) do taskid, abort
                if taskid == 1
                    error("boom from 1")
                else
                    # Deterministic: wait on the very signal we need, not a guessed delay.
                    ok = wait_for(() -> abort[])
                    abort_observed[taskid - 1] = ok
                    Threads.atomic_add!(done, 1)
                end
            end
        catch caught
            e = caught
        end
        @test e isa TaskFailedException
        @test e.task.exception isa ErrorException
        @test e.task.exception.msg == "boom from 1"
        @test all(abort_observed) # every sibling actually saw `abort[] == true`
        @test done[] == 3 # ... and every sibling was actually waited for before returning
    end

    @testset "first *actual* failure is reported, not the lowest taskid" begin
        # taskid 3 is made to fail before taskid 1, despite 1 < 3. Task 1 waits for `abort[]`
        # itself (not a separate flag) before failing: `abort[]` is only ever set *inside* a
        # catch block, right next to the CAS that records the failure, so by the time task 1
        # observes it, task 3's own catch block (including its CAS) has already run - no
        # sleep-based margin needed, this is the actual guarantee being tested, used as its own
        # synchronization. This directly tests that `spawn_tasks_safe_wait` records failures by
        # actual order, not by iterating `1:num_tasks`.
        e = nothing
        try
            FerriteAssembly.spawn_tasks_safe_wait(4) do taskid, abort
                if taskid == 3
                    error("boom from 3")
                elseif taskid == 1
                    wait_for(() -> abort[]) || error("handshake timed out")
                    error("boom from 1")
                else
                    wait_for(() -> abort[])
                end
            end
        catch caught
            e = caught
        end
        @test e isa TaskFailedException
        @test e.task.exception isa ErrorException
        @test e.task.exception.msg == "boom from 3"
    end

    @testset "a child's own InterruptException is reported like any other failure" begin
        # Distinguishes a task failing *with* an InterruptException (just another error, as
        # far as spawn_tasks_safe_wait is concerned) from an InterruptException delivered to the
        # *caller* while it waits (tested below) - the latter must not be confused with the
        # former by the drain loop's `catch e; e isa InterruptException` check, since `wait`
        # on a task that failed throws `TaskFailedException`, not the bare `InterruptException`.
        e = nothing
        try
            FerriteAssembly.spawn_tasks_safe_wait(2) do taskid, abort
                if taskid == 1
                    throw(InterruptException())
                else
                    wait_for(() -> abort[])
                end
            end
        catch caught
            e = caught
        end
        @test e isa TaskFailedException
        @test e.task.exception isa InterruptException
    end

    @testset "InterruptException delivered to the waiting caller does not skip draining" begin
        # Task 1 cannot finish until `release[]` is set - its completion is gated explicitly,
        # not by a guessed sleep duration, so there is no way for it to race ahead of the
        # interrupt. The *interruptor* task both delivers the interrupt and (after a short,
        # one-directional grace period - not a guess about *this* test's correctness, only
        # about how long an already-sent, same-process async exception takes to land) sets
        # `release[]` itself, in that order on a single task, so "interrupt sent" strictly
        # precedes "task 1 allowed to finish" with no race between those two steps either.
        started = Threads.Atomic{Bool}(false)
        release = Threads.Atomic{Bool}(false)
        finished = Threads.Atomic{Bool}(false)
        maintask = current_task()
        interruptor = Threads.@spawn begin
            if wait_for(() -> started[])
                schedule(maintask, InterruptException(); error = true)
                sleep(0.05) # let the already-sent exception actually land in `maintask`
            end
            release[] = true
        end
        caught = nothing
        try
            FerriteAssembly.spawn_tasks_safe_wait(2) do taskid, abort
                if taskid == 1
                    started[] = true
                    wait_for(() -> release[]; timeout = 30.0)
                    finished[] = true
                end
            end
        catch e
            caught = e
        end
        # Snapshot `finished[]` the instant `spawn_tasks_safe_wait` returns/throws, *before* waiting
        # for the interruptor below: that wait only returns once `release[]` is set, which lets
        # task 1 finish independently of whether `spawn_tasks_safe_wait` itself waited correctly - so
        # reading `finished[]` afterwards would pass even for a broken helper that returns as
        # soon as the interrupt arrives (confirmed: a deliberately-broken in-memory variant,
        # with the interrupt retry's `continue` replaced by `break`, passed this assertion in
        # 3 of 5 runs when checked only after that wait).
        finished_at_return = finished[]
        # Defensive: a timed-out handshake above would make the interruptor task fail (not
        # throw into `maintask`); don't let that crash the whole test suite.
        try
            wait(interruptor)
        catch
        end
        @test caught isa InterruptException
        @test finished_at_return # the slow task had already finished by the time the call returned
    end

    @testset "an unrelated TaskFailedException delivered asynchronously does not skip draining" begin
        # A `TaskFailedException` isn't only produced by `wait(t)` discovering that `t` itself
        # failed - it can also arrive *asynchronously* (e.g. forwarded by external supervising
        # code watching some other task), exactly like the `InterruptException` case above. The
        # drain loop must not mistake "I caught a TaskFailedException" for "the task I'm
        # currently waiting on has failed": it must check which task the exception actually
        # refers to, or it stops waiting on a sibling that is still running.
        other = Threads.@spawn error("unrelated failure, not one of spawn_tasks_safe_wait's own tasks")
        try
            wait(other)
        catch
        end # let `other` actually fail before using it below

        started = Threads.Atomic{Bool}(false)
        release = Threads.Atomic{Bool}(false)
        finished = Threads.Atomic{Bool}(false)
        maintask = current_task()
        interruptor = Threads.@spawn begin
            if wait_for(() -> started[])
                schedule(maintask, TaskFailedException(other); error = true)
                sleep(0.05) # let the already-sent exception actually land in `maintask`
            end
            release[] = true
        end
        caught = nothing
        try
            FerriteAssembly.spawn_tasks_safe_wait(2) do taskid, abort
                if taskid == 1
                    started[] = true
                    wait_for(() -> release[]; timeout = 30.0)
                    finished[] = true
                end
            end
        catch e
            caught = e
        end
        # Snapshot immediately, before waiting on `interruptor` below, for the same reason as
        # the `InterruptException` test above: that wait only returns once `release[]` is set,
        # which would let task 1 finish independently of whether draining was done correctly.
        finished_at_return = finished[]
        try
            wait(interruptor)
        catch
        end
        @test caught isa TaskFailedException
        @test caught.task === other
        @test finished_at_return # the slow task had already finished by the time the call returned
    end
end

# `drain_tasks` is the lower-level primitive `spawn_tasks_safe_wait` builds on; it is only ever
# called with freshly `Threads.@spawn`ed tasks in practice, never the calling task itself, so
# this case cannot arise via `spawn_tasks_safe_wait`'s public behavior. Tested directly here
# since there is no way to reach it through the public API.
@testset "drain_tasks rethrows immediately (not an infinite retry) for a self-wait" begin
    # `wait(t; throw = false)`'s self-wait guard (`t === current_task()`) is stateless and
    # unconditional - it throws `ConcurrencyViolationError` the same way on every call, so
    # treating it like a retryable asynchronous exception would spin forever rather than ever
    # returning. If this test hangs instead of completing, that regression is exactly what it
    # is meant to catch.
    @test_throws ConcurrencyViolationError FerriteAssembly.drain_tasks([current_task()])
end
