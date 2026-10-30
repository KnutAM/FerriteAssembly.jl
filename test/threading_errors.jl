# Regression tests for the race condition where `work_domain_threaded!` (src/work.jl) used
# `Base.Experimental.@sync`, which rethrows the first failing task's exception *without*
# waiting for the other spawned tasks. Those tasks kept running on the same persistent
# task-local item buffers (and kept writing into the same `K`/`r`) as a subsequent `work!` call,
# e.g. a caller that catches the exception, halves the time step, and reassembles.
#
# The hard concurrency guarantees themselves (no orphaned tasks, true first-failure order,
# interrupt-safety) are tested directly and fast against the extracted primitive,
# `FerriteAssembly.sync_with_abort`, in test/task_synchronization.jl - that needs no grid or
# CellValues, only trivial closures, so it isn't sensitive to scheduling noise under load.
# This file only checks that `work_domain_threaded!` actually *uses* that primitive correctly:
# a failing element routine surfaces as `TaskFailedException`, assembly stops early rather than
# continuing through the whole grid, and a retry afterwards (reusing the same persistent
# task-local buffers) matches a sequential reference - both directly and through
# `CoupledSimulations`.
@testset "Threaded work! finishes all tasks before propagating a failure (integration)" begin
    num_tasks = 4 # deliberately independent of `JULIA_NUM_THREADS` (CI runs with 2)
    # Many cells/chunks, so that however many extra chunks siblings race through before
    # noticing `abort[]` (there is no artificial delay here, unlike the old barrier-based
    # version of this test), "early abort" vs. "processed (almost) everything" stays an
    # unmistakable, non-flaky distinction - see the assertion below.
    grid = generate_grid(Quadrilateral, (50, 50))
    ip = Lagrange{RefQuadrilateral, 1}()
    dh = close!(add!(DofHandler(grid), :u, ip))
    cv = CellValues(QuadratureRule{RefQuadrilateral}(2), ip)
    throwcell = 1

    FAILING = Threads.Atomic{Bool}(false) # true during the call that should fail
    CALLS = Threads.Atomic{Int}(0)        # total element routine invocations

    struct RaceMaterial end

    function FerriteAssembly.element_routine!(Ke, re, state, ae, ::RaceMaterial, cv, buffer)
        Threads.atomic_add!(CALLS, 1)
        if FAILING[] && cellid(buffer) == throwcell
            error("injected failure in cell $throwcell")
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
        return nothing
    end

    db = setup_domainbuffer(DomainSpec(dh, RaceMaterial(), cv); threading = true, num_tasks)
    db_ref = setup_domainbuffer(DomainSpec(dh, RaceMaterial(), cv); threading = false)

    a = zeros(ndofs(dh))
    K_ref = allocate_matrix(dh)
    r_ref = zeros(ndofs(dh))
    work!(KeReAssembler(K_ref, r_ref), db_ref; a)

    @testset "failure surfaces, assembly stops early, retry matches reference" begin
        FAILING[] = true
        calls_before = CALLS[]
        K = allocate_matrix(dh)
        r = zeros(ndofs(dh))
        @test_throws TaskFailedException work!(KeReAssembler(K, r), db; a)
        # `sync_with_abort` has already fully drained every task by the time `work!` returns
        # (or throws), so this is stable right away - no sleeping/polling needed.
        @test CALLS[] - calls_before < getncells(grid) ÷ 2

        FAILING[] = false
        work!(KeReAssembler(K, r), db; a) # retry: fresh KeReAssembler re-zeros K, r
        @test K ≈ K_ref
        @test r ≈ r_ref
    end

    @testset "through CoupledSimulations" begin
        dh2 = close!(add!(DofHandler(grid), :v, ip))
        struct RaceTrivialPartner end
        FerriteAssembly.element_routine!(Ke, re, state, ae, ::RaceTrivialPartner, cv, buffer) = nothing
        db2 = setup_domainbuffer(DomainSpec(dh2, RaceTrivialPartner(), cv); threading = true, num_tasks)
        sim1 = Simulation(db, a)
        sim2 = Simulation(db2, zeros(ndofs(dh2)))
        g = CoupledSimulations((u = sim1,); refs = (v = sim2,))

        FAILING[] = true
        calls_before = CALLS[]
        K = allocate_matrix(dh)
        r = zeros(ndofs(dh))
        @test_throws TaskFailedException work!(KeReAssembler(K, r), g.u)
        @test CALLS[] - calls_before < getncells(grid) ÷ 2

        FAILING[] = false
        work!(KeReAssembler(K, r), g.u) # retry must succeed and match the sequential reference
        @test K ≈ K_ref
        @test r ≈ r_ref
    end
end
