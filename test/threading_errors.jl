# Regression tests for the race condition where `work_domain_threaded!` (src/work.jl) used
# `Base.Experimental.@sync`, which rethrows the first failing task's exception *without*
# waiting for the other spawned tasks. Those tasks kept running on the same persistent
# task-local item buffers (and kept writing into the same `K`/`r`) as a subsequent `work!` call,
# e.g. a caller that catches the exception, halves the time step, and reassembles.
#
# The hard concurrency guarantees themselves (no orphaned tasks, true first-failure order,
# interrupt-safety) are tested directly and fast against the extracted primitive,
# `FerriteAssembly.spawn_tasks_safe_wait`, in test/task_synchronization.jl - that needs no grid or
# CellValues, only trivial closures, so it isn't sensitive to scheduling noise under load.
# This file only checks that `work_domain_threaded!` actually *uses* that primitive correctly:
# a failing element routine surfaces as `TaskFailedException`, assembly stops early rather than
# continuing through the whole grid, and a retry afterwards (reusing the same persistent
# task-local buffers) matches a sequential reference - both directly and through
# `CoupledSimulations`.
@testset "Threaded work! finishes all tasks before propagating a failure (integration)" begin
    num_tasks = 4 # deliberately independent of `JULIA_NUM_THREADS` (CI runs with 2)
    grid = generate_grid(Quadrilateral, (50, 50))
    ip = Lagrange{RefQuadrilateral, 1}()
    dh = close!(add!(DofHandler(grid), :u, ip))
    cv = CellValues(QuadratureRule{RefQuadrilateral}(2), ip)
    # `work_domain_threaded!` processes one color (one `spawn_tasks_safe_wait` call) at a time, and
    # a color earlier than the throwing one is *always* fully processed regardless of whether
    # abort checks work (the outer loop over colors simply never reaches later colors once one
    # call throws) - so bounding "cells processed after the failing call" against the *whole*
    # grid would pass even with the abort checks in the assembly loop removed entirely, as long
    # as the throwing color plus every earlier one summed to under half the grid. Picking the
    # throw cell from the *first* color and bounding against just that color's size avoids that
    # loophole: with the checks removed, every cell in this color would run (not a handful).
    colors = create_coloring(grid)
    # `setup_domainbuffer` sorts each color before splitting it into chunks (see
    # `create_chunks`/`split_in_chunks` in src/Multithreading/TaskChunks.jl), so the cell that
    # actually ends up first in the first chunk is the *smallest* cell number in `colors[1]`,
    # not whatever `create_coloring` happens to return first (unsorted).
    throwcell = minimum(colors[1])

    FAILING = Threads.Atomic{Bool}(false) # true during the call that should fail
    CALLS = Threads.Atomic{Int}(0)        # total element routine invocations

    struct RaceMaterial end

    function FerriteAssembly.element_routine!(Ke, re, state, ae, ::RaceMaterial, cv, buffer)
        Threads.atomic_add!(CALLS, 1)
        if FAILING[]
            if cellid(buffer) == throwcell
                error("injected failure in cell $throwcell")
            end
            # Under true multi-threading there is nothing else synchronizing "the throw"
            # against "how much a sibling can race through in the meantime" - with a trivial,
            # no-delay element routine, a sibling can process most or all of the remaining
            # color before ever checking `abort[]` again, independent of whether the abort
            # checks in the assembly loop are even present (confirmed empirically: this
            # happened in a majority of trials without this delay). A tiny fixed delay, applied
            # only during the (two) failing calls in this file, makes "how much extra work can
            # happen before `abort[]` is noticed" bounded by that delay rather than by
            # unpredictable scheduling - no computation needed here since K/r correctness is
            # only checked via the separate, non-failing retry call below.
            sleep(0.001)
            return nothing
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
    # Warm up (compile) the threaded dispatch chain before the failing round below: compiling
    # it for the first time *during* that round (a fresh process, as in CI, has not compiled
    # any of this yet) would give siblings extra wall-clock time to race through the throwing
    # color before noticing `abort[]`, undermining the early-abort assertion.
    work!(KeReAssembler(allocate_matrix(dh), zeros(ndofs(dh))), db; a)

    @testset "failure surfaces, assembly stops early, retry matches reference" begin
        FAILING[] = true
        calls_before = CALLS[]
        K = allocate_matrix(dh)
        r = zeros(ndofs(dh))
        @test_throws TaskFailedException work!(KeReAssembler(K, r), db; a)
        # `spawn_tasks_safe_wait` has already fully drained every task by the time `work!` returns
        # (or throws), so this is stable right away - no sleeping/polling needed.
        @test CALLS[] - calls_before < length(colors[1]) ÷ 2

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
        # Warm up (compile) the coupled dispatch chain (reinit_buffer!/work_single! through
        # CoupledSimulation) before the failing round: compiling it for the first time *during*
        # that round would give siblings extra wall-clock time to race through the throwing
        # color before noticing `abort[]`, undermining the early-abort assertion below.
        work!(KeReAssembler(allocate_matrix(dh), zeros(ndofs(dh))), g.u)

        FAILING[] = true
        calls_before = CALLS[]
        K = allocate_matrix(dh)
        r = zeros(ndofs(dh))
        @test_throws TaskFailedException work!(KeReAssembler(K, r), g.u)
        @test CALLS[] - calls_before < length(colors[1]) ÷ 2

        FAILING[] = false
        work!(KeReAssembler(K, r), g.u) # retry must succeed and match the sequential reference
        @test K ≈ K_ref
        @test r ≈ r_ref
    end
end
