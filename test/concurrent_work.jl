# Tests for the question: with PR110's fix (work_domain_threaded! always finishes every
# spawned task before propagating a failure or returning), is it safe for independent
# *outer* callers - e.g. separate points of a parameter sweep run under `Threads.@spawn` -
# to call threaded `work!` concurrently, each with its own inner `spawn_tasks_safe_wait`
# tasks (genuine two-level nested threading)?
#
# Answer exercised here: yes, as long as each outer call works its own, independently
# `setup_domainbuffer`'d buffer (own K/r, own states, own item-buffer scratch). See
# test/replacements.jl's "replace_material shares item-buffer scratch..." testset for the
# one common way two buffers end up *not* independent despite looking like separate objects.
#
# Overlap between the two outer calls is confirmed, not assumed: each element-routine
# invocation records its wall-clock time, and a small fixed per-item delay elongates each
# outer call enough that the two recorded [first call, last call] intervals can be checked
# to genuinely overlap. This avoids a live, blocking handshake between the two outer calls
# (tried first; abandoned because busy-waiting inside the hot assembly loop under this much
# task oversubscription could itself stall for many seconds) while still *proving* the two
# outer calls were truly running at the same time rather than merely both being `@spawn`ed.
@testset "Concurrent threaded work! on independent domain buffers (parameter sweep)" begin
    grid = generate_grid(Quadrilateral, (20, 20))
    ip = Lagrange{RefQuadrilateral,1}()
    dh = close!(add!(DofHandler(grid), :u, ip))
    cv = CellValues(QuadratureRule{RefQuadrilateral}(2), ip)
    ks = (1.0, 2.0) # the two "parameter values" being swept concurrently
    num_tasks = 2    # each outer call's own inner task count, independent of Threads.nthreads()

    t_start = (Ref(0.0), Ref(0.0))
    t_end = (Ref(0.0), Ref(0.0))
    timing_lock = ReentrantLock()
    fail_on = (Threads.Atomic{Int}(0), Threads.Atomic{Int}(0))

    struct SweepTag{N} end

    function FerriteAssembly.element_routine!(Ke, re, state, ae, ::SweepTag{N}, cv, buffer) where N
        sleep(1e-4) # elongates each outer call enough to reliably overlap its sibling, see above
        now = time()
        lock(timing_lock) do
            t_start[N][] == 0.0 && (t_start[N][] = now)
            t_end[N][] = now
        end
        fc = fail_on[N][]
        if fc != 0 && cellid(buffer) == fc
            # Bound how far a sibling *inner* task of this same outer call can race ahead
            # before observing abort[], see test/threading_errors.jl for the same technique.
            sleep(0.001)
            error("injected failure in sweep $N on cell $fc")
        end
        FerriteAssembly.element_routine!(Ke, re, state, ae, EE.StationaryFourier(ks[N]), cv, buffer)
    end

    dbs = ntuple(N -> setup_domainbuffer(DomainSpec(dh, SweepTag{N}(), cv); threading=true, num_tasks), 2)
    dbs_ref = ntuple(N -> setup_domainbuffer(DomainSpec(dh, EE.StationaryFourier(ks[N]), cv); threading=false), 2)

    # Derived from `dbs[1]`'s own chunking, rather than independently recomputing
    # `create_coloring(grid)` in this file: the two are not guaranteed to agree on which
    # cells land in "color 1", since greedy coloring depends on the iteration order of the
    # cellset it is given, and `setup_domainbuffer` colors a `Set`, not a `UnitRange`.
    # `chunks[1]` is `dbs[1]`'s actual first color, split into per-task chunks; cells within
    # it are sorted, so its smallest cell number is the first one any task reaches.
    color1_cells = reduce(vcat, first(FerriteAssembly.get_chunks(dbs[1])))
    throwcell = minimum(color1_cells)

    # The two sweep buffers are genuinely independent: no shared item-buffer scratch or
    # state, unlike the `replace_material` case covered in test/replacements.jl.
    ib1, ib2 = FerriteAssembly.get_itembuffer(dbs[1]), FerriteAssembly.get_itembuffer(dbs[2])
    @test FerriteAssembly.get_local(ib1, 1).Ke !== FerriteAssembly.get_local(ib2, 1).Ke
    @test FerriteAssembly.get_state(dbs[1]) !== FerriteAssembly.get_state(dbs[2])

    a = zeros(ndofs(dh))
    Ks_ref = (allocate_matrix(dh), allocate_matrix(dh))
    rs_ref = (zeros(ndofs(dh)), zeros(ndofs(dh)))
    for N in 1:2
        work!(KeReAssembler(Ks_ref[N], rs_ref[N]), dbs_ref[N]; a)
    end
    # Warm up (compile) the dispatch chain before measuring, as in test/threading_errors.jl.
    for N in 1:2
        work!(KeReAssembler(allocate_matrix(dh), zeros(ndofs(dh))), dbs[N]; a)
    end

    @testset "both complete correctly, with confirmed overlap" begin
        for N in 1:2
            t_start[N][] = 0.0
            t_end[N][] = 0.0
        end
        Ks = (allocate_matrix(dh), allocate_matrix(dh))
        rs = (zeros(ndofs(dh)), zeros(ndofs(dh)))
        tasks = [Threads.@spawn(work!(KeReAssembler(Ks[N], rs[N]), dbs[N]; a)) for N in 1:2]
        foreach(wait, tasks) # propagates any unexpected exception immediately

        # The two outer `work!` calls' recorded [first call, last call] wall-clock
        # intervals genuinely overlap - not just both having been `@spawn`ed.
        @test t_start[1][] < t_end[2][] && t_start[2][] < t_end[1][]
        for N in 1:2
            @test Ks[N] ≈ Ks_ref[N]
            @test rs[N] ≈ rs_ref[N]
        end
    end

    @testset "a failure in one sweep point doesn't affect the concurrently running sibling" begin
        # Genuine overlap between the two outer calls was already established above, using
        # the same buffers and the same per-item delay; it is not re-checked here against a
        # wall-clock bound, since injecting the failure means sweep 1 can legitimately finish
        # in very few calls, making a timing-based overlap (or early-stop count) assertion
        # flaky under an unlucky scheduling order rather than a meaningful correctness check.
        # The "stops early rather than finishing the whole grid" guarantee itself is already
        # covered, without any such timing sensitivity, by the dedicated primitive tests in
        # test/task_synchronization.jl and by the single-outer-call integration test in
        # test/threading_errors.jl. What's new and worth asserting here is specific to the
        # *nested* case: the failure's shape/identity, and that it leaves the concurrently
        # running sibling outer call entirely correct and unaffected.
        fail_on[1][] = throwcell
        K1, r1 = allocate_matrix(dh), zeros(ndofs(dh))
        K2, r2 = allocate_matrix(dh), zeros(ndofs(dh))
        t1 = Threads.@spawn work!(KeReAssembler(K1, r1), dbs[1]; a)
        t2 = Threads.@spawn work!(KeReAssembler(K2, r2), dbs[2]; a)
        e1 = nothing
        try
            wait(t1)
        catch caught
            e1 = caught
        end
        wait(t2) # must complete successfully and correctly despite the sibling's failure
        fail_on[1][] = 0 # reset, in case this testset is ever re-run within the same session

        # `work!` itself throws `TaskFailedException(inner_task)` (from `spawn_tasks_safe_wait`),
        # so the outer `@spawn`ed task's own stored exception is a `TaskFailedException`
        # wrapping the actually-failing *inner* task, not the `ErrorException` directly.
        @test e1 isa TaskFailedException
        @test e1.task.exception isa TaskFailedException
        inner_exc = e1.task.exception.task.exception
        @test inner_exc isa ErrorException
        @test inner_exc.msg == "injected failure in sweep 1 on cell $throwcell"
        @test K2 ≈ Ks_ref[2]
        @test r2 ≈ rs_ref[2]
    end
end
