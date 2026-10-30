# Regression tests for the race condition where `work_domain_threaded!` (src/work.jl) used
# `Base.Experimental.@sync`, which rethrows the first failing task's exception *without*
# waiting for the other spawned tasks. Those tasks kept running on the same persistent
# task-local item buffers (and kept writing into the same `K`/`r`) as a subsequent `work!` call,
# e.g. a caller that catches the exception, halves the time step, and reassembles.
#
# The fault is injected deterministically: with `num_tasks` tasks and explicit `chunks` giving
# the first color at least `num_tasks` chunks, every task's first item blocks (one of them is
# the designated cell that throws, after waiting for the others to reach the same barrier). The
# others stay blocked until released by a bounded timer, so items already in progress get a
# chance to finish before a fixed `work!` returns - mirroring the upstream reproducer this is
# based on.
@testset "Threaded work! finishes all tasks before propagating a failure" begin
    num_tasks = 4 # deliberately independent of `JULIA_NUM_THREADS` (CI runs with 2)
    grid = generate_grid(Quadrilateral, (10, 10))
    ip = Lagrange{RefQuadrilateral, 1}()
    dh = close!(add!(DofHandler(grid), :u, ip))
    cv = CellValues(QuadratureRule{RefQuadrilateral}(2), ip)

    colors = create_coloring(grid)
    chunk_size = 2
    chunks = [collect(Iterators.partition(sort(c), chunk_size)) .|> collect for c in colors]
    @test length(first(chunks)) >= num_tasks # every task must get a chunk of the first color
    throwcell = first(first(first(chunks)))
    nbarrier = min(num_tasks, length(first(chunks)))

    FAILING = Threads.Atomic{Bool}(false)  # true during the call that should fail
    RELEASE = Threads.Atomic{Bool}(false)  # releases the tasks blocked in the failing call
    ENTERED = Threads.Atomic{Int}(0)       # non-throwing tasks that reached the barrier
    ACTIVE = Threads.Atomic{Int}(0)        # element routines currently executing
    CALLS = Threads.Atomic{Int}(0)         # total element routine invocations
    BARRIER_OK = Threads.Atomic{Bool}(true) # false if any `wait_for` below timed out
    TIMER = Ref{Union{Timer, Nothing}}(nothing)
    SLOW = Threads.Atomic{Bool}(false)      # true: `throwcell` sleeps instead of failing
    SLOWDONE = Threads.Atomic{Bool}(false)  # set once that sleep has completed

    function wait_for(cond; timeout = 30.0)
        t0 = time()
        while !cond()
            time() - t0 > timeout && return false
            # A bare `yield()` only reschedules already-queued tasks; it does not reliably let
            # the libuv event loop (which drives `RELEASE`'s `Timer` callback below) run when
            # every `:default`-pool thread is busy in this same loop. `sleep` does (it is built
            # on a `Timer` + `wait` itself), so poll with a short sleep instead.
            sleep(0.001)
        end
        return true
    end

    struct RaceMaterial end

    function FerriteAssembly.element_routine!(Ke, re, state, ae, ::RaceMaterial, cv, buffer)
        Threads.atomic_add!(CALLS, 1)
        Threads.atomic_add!(ACTIVE, 1)
        try
            if SLOW[] && cellid(buffer) == throwcell
                sleep(0.5)
                SLOWDONE[] = true
            elseif FAILING[]
                if cellid(buffer) == throwcell
                    wait_for(() -> ENTERED[] >= nbarrier - 1) || (BARRIER_OK[] = false)
                    # Release the others even if no retry starts (needed to let a *fixed*
                    # `work!` call return instead of blocking forever on its own siblings).
                    # A `Base.Timer` (libuv-backed, like `sleep`) is used rather than another
                    # `Threads.@spawn`ed task: with `num_tasks == Threads.nthreads()`, all
                    # `:default`-pool threads are already pinned running the busy-waiting
                    # siblings below, so a competing spawned task could be starved indefinitely
                    # waiting for a free thread - a `Timer`'s callback does not need one.
                    TIMER[] = Timer(_ -> (RELEASE[] = true), 0.3)
                    error("injected failure in cell $throwcell")
                elseif !RELEASE[]
                    Threads.atomic_add!(ENTERED, 1)
                    wait_for(() -> RELEASE[]) || (BARRIER_OK[] = false)
                end
            end
            for q in 1:getnquadpoints(cv)
                dΩ = getdetJdV(cv, q)
                for i in 1:getnbasefunctions(cv)
                    δN = shape_value(cv, q, i)
                    ∇δN = shape_gradient(cv, q, i)
                    re[i] -= δN * dΩ
                    for j in 1:getnbasefunctions(cv)
                        Ke[i, j] += (∇δN ⋅ shape_gradient(cv, q, j)) * dΩ
                    end
                end
            end
            re .+= Ke * ae
        finally
            Threads.atomic_sub!(ACTIVE, 1)
        end
        return nothing
    end

    db = setup_domainbuffer(DomainSpec(dh, RaceMaterial(), cv; chunks); threading = true, num_tasks)
    db_ref = setup_domainbuffer(DomainSpec(dh, RaceMaterial(), cv); threading = false)

    a = zeros(ndofs(dh))
    K_ref = allocate_matrix(dh)
    r_ref = zeros(ndofs(dh))
    work!(KeReAssembler(K_ref, r_ref), db_ref; a)
    # Warm up (compile) the full threaded dispatch chain for `RaceMaterial` outside of the
    # fault-injection rounds below: the barrier logic assumes `ENTERED`/`RELEASE` are observed
    # promptly, which first-time JIT compilation (done concurrently by several tasks the first
    # time `work!` threads through a new material) could otherwise delay unpredictably.
    work!(KeReAssembler(allocate_matrix(dh), zeros(ndofs(dh))), db; a)

    function reset_round!()
        FAILING[] = false; RELEASE[] = false; ENTERED[] = 0; BARRIER_OK[] = true; TIMER[] = nothing
        SLOW[] = false; SLOWDONE[] = false
    end
    function drain!()
        t = TIMER[]
        # `wait` on an already-fired one-shot `Timer` throws `EOFError` (unlike `Task`'s
        # `wait`, which just returns); `isopen` tells us whether it has fired yet.
        t === nothing || !isopen(t) || wait(t)
        sleep(0.5) # let any still-running (pre-fix: orphaned) tasks finish
    end

    @testset "no orphaned tasks, and work! stops making progress" begin
        reset_round!()
        FAILING[] = true
        calls_before = CALLS[]
        K = allocate_matrix(dh)
        r = zeros(ndofs(dh))
        @test_throws TaskFailedException work!(KeReAssembler(K, r), db; a)
        # Right when `work!` has thrown, every spawned task (the failing one and the ones it
        # was racing with) must already be done - none left running in the background.
        @test ACTIVE[] == 0
        @test BARRIER_OK[]
        calls_at_throw = CALLS[]
        @test calls_at_throw - calls_before == num_tasks # throwing cell + released siblings
        sleep(0.5) # longer than the release timer
        # No element routine may start after `work!` has returned: early abort, not just
        # "eventually catches up" - far fewer than the grid's `getncells(grid)` cells.
        @test CALLS[] == calls_at_throw
        drain!()
    end

    @testset "retry after a failure matches the sequential reference" begin
        reset_round!()
        FAILING[] = true
        K = allocate_matrix(dh)
        r = zeros(ndofs(dh))
        @test_throws TaskFailedException work!(KeReAssembler(K, r), db; a)
        drain!()
        reset_round!()
        work!(KeReAssembler(K, r), db; a) # fresh KeReAssembler re-zeros K, r
        @test K ≈ K_ref
        @test r ≈ r_ref
    end

    @testset "exception type and content are unchanged" begin
        reset_round!()
        FAILING[] = true
        K = allocate_matrix(dh)
        r = zeros(ndofs(dh))
        e = nothing
        try
            work!(KeReAssembler(K, r), db; a)
        catch caught
            e = caught
        end
        @test e isa TaskFailedException
        @test e.task.exception isa ErrorException
        @test occursin("injected failure in cell $throwcell", e.task.exception.msg)
        drain!()
    end

    @testset "an InterruptException while waiting does not skip draining" begin
        # Regression test for a narrower variant of the same bug: `work_domain_threaded!`'s
        # drain loop must not treat an `InterruptException` delivered to the *waiting* caller
        # (e.g. Ctrl-C) as if the task being waited on had finished - that task may still be
        # running, and `work!` must still not return until it actually has.
        reset_round!()
        SLOW[] = true
        maintask = current_task()
        interruptor = Threads.@spawn begin
            sleep(0.1) # give `work!` time to start and block in `wait` on the slow task
            schedule(maintask, InterruptException(); error = true)
        end
        K = allocate_matrix(dh)
        r = zeros(ndofs(dh))
        caught = nothing
        try
            work!(KeReAssembler(K, r), db; a)
        catch e
            caught = e
        end
        wait(interruptor)
        @test caught isa InterruptException
        @test SLOWDONE[] # the slow task had already finished by the time `work!` returned
        reset_round!()
    end

    @testset "through CoupledSimulations" begin
        reset_round!()
        dh2 = close!(add!(DofHandler(grid), :v, ip))
        struct RaceTrivialPartner end
        FerriteAssembly.element_routine!(Ke, re, state, ae, ::RaceTrivialPartner, cv, buffer) = nothing
        db2 = setup_domainbuffer(DomainSpec(dh2, RaceTrivialPartner(), cv); threading = true, num_tasks)
        sim1 = Simulation(db, a)
        sim2 = Simulation(db2, zeros(ndofs(dh2)))
        g = CoupledSimulations((u = sim1,); refs = (v = sim2,))
        work!(KeReAssembler(allocate_matrix(dh), zeros(ndofs(dh))), g.u) # warm up (compile)

        FAILING[] = true
        calls_before = CALLS[]
        K = allocate_matrix(dh)
        r = zeros(ndofs(dh))
        @test_throws TaskFailedException work!(KeReAssembler(K, r), g.u)
        @test ACTIVE[] == 0
        @test BARRIER_OK[]
        calls_at_throw = CALLS[]
        sleep(0.5)
        @test CALLS[] == calls_at_throw # no orphaned task kept going through the coupled path
        drain!()

        reset_round!()
        work!(KeReAssembler(K, r), g.u) # retry must succeed and match the sequential reference
        @test K ≈ K_ref
        @test r ≈ r_ref
    end
end
