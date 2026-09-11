@testset "replace_material" begin
    m_el = EE.LinearElastic(;E=1.0, ν=0.4)
    m_elx2 = EE.LinearElastic(;E=2.0, ν=0.4)
    f_repl1(::EE.LinearElastic) = m_elx2
    m_pl = EE.J2Plasticity(;E=1.0, ν=0.4, σ0=0.2, H=1.0)
    f_repl2(::EE.LinearElastic) = m_pl
    f_repl3(::EE.J2Plasticity) = m_el
    grid = generate_grid(Quadrilateral, (2,2))
    ip = Lagrange{RefQuadrilateral,1}()^2
    dh = DofHandler(grid); add!(dh, :u, ip); close!(dh)
    qr = QuadratureRule{RefQuadrilateral}(2)
    cv = CellValues(qr, ip, ip)
    dspec = DomainSpec(dh, m_el, cv)
    buffer = setup_domainbuffer(dspec)
    ad_buffer = setup_domainbuffer(dspec; autodiffbuffer=true)
    td_buffer = setup_domainbuffer(dspec; threading=true)
    
    for b0 in (buffer, ad_buffer, td_buffer)
        @test FerriteAssembly.get_material(b0) === m_el
        b1 = FerriteAssembly.replace_material(b0, f_repl1)
        @test FerriteAssembly.get_material(b1) === m_elx2
        b2 = FerriteAssembly.replace_material(b1, f_repl2)
        @test FerriteAssembly.get_material(b2) === m_pl
        b3 = FerriteAssembly.replace_material(b2, f_repl3)
        @test FerriteAssembly.get_material(b3) === m_el
    end

    n_half = getncells(grid)÷2
    buffers = setup_domainbuffers(Dict(
        "a" => DomainSpec(dh, m_el, cv; set=1:(n_half-1)),
        "b" => DomainSpec(dh, m_pl, cv; set=n_half:getncells(grid))
    ))
    @test FerriteAssembly.get_material(buffers, "a") === m_el
    @test FerriteAssembly.get_material(buffers, "b") === m_pl 
    f_repl(::EE.LinearElastic) = m_elx2
    f_repl(m::EE.J2Plasticity) = m
    bs2 = FerriteAssembly.replace_material(buffers, f_repl)
    @test FerriteAssembly.get_material(bs2, "a") === m_elx2
    @test FerriteAssembly.get_material(bs2, "b") === m_pl
end

@testset "couple_buffers" begin
    grid = generate_grid(Quadrilateral, (2,2))
    addcellset!(grid, "left", x -> x[1] < eps())
    addcellset!(grid, "right", setdiff(1:getncells(grid), getcellset(grid, "left")))
    ip = Lagrange{RefQuadrilateral,1}()
    dh1 = close!(add!(DofHandler(grid), :u, ip))
    dh2 = close!(add!(DofHandler(grid), :v, ip^2))
    qr = QuadratureRule{RefQuadrilateral}(2)
    cvu = CellValues(qr, ip, ip)
    cvv = CellValues(qr, ip^2, ip)

    struct MA end
    struct MB end
    # We will test with the following dof value differences 
    # aold will be same in both cases (for both components in the case of MB)
    # a will be 3 times larger for first component in MB, and 5 times for second component
    # State will be 6 times larger for MB, obtained by multiplying the function values by factor 2
    FerriteAssembly.create_cell_state(::MA, cv, x, ae, args...) = [function_value(cv, i, ae) for i in 1:getnquadpoints(cv)]
    FerriteAssembly.create_cell_state(::MB, cv, x, ae, args...) = [2 * function_value(cv, i, ae)[1] for i in 1:getnquadpoints(cv)]

    Δt2 = 0.25
    # Test case to check that values have been updated correctly
    function FerriteAssembly.element_routine!(Ke, re, state, ae, m::MA, cv, buffer)
        cb_b = FerriteAssembly.get_coupled_buffer(buffer, :b)
        # Check that correct material has been set
        @test FerriteAssembly.get_material(cb_b) isa MB
        # Check that dofs have been updated
        @test 3 * ae ≈ FerriteAssembly.get_ae(cb_b)[1:2:end] # 1st component
        @test 5 * ae ≈ FerriteAssembly.get_ae(cb_b)[2:2:end] # 2nd component
        # Check that old dofs have been updated
        @test FerriteAssembly.get_aeold(buffer) ≈ FerriteAssembly.get_aeold(cb_b)[1:2:end]
        @test FerriteAssembly.get_aeold(buffer) ≈ FerriteAssembly.get_aeold(cb_b)[2:2:end]
        # Check that state variables have been updated
        @test 6 * state ≈ FerriteAssembly.get_state(cb_b)
        # Check that the coupled buffer's time increment reflects the partner's current value
        @test FerriteAssembly.get_time_increment(cb_b) == Δt2
    end

    a1 = rand(ndofs(dh1))
    a2 = zeros(ndofs(dh2))
    @assert length(a1) * 2 == length(a2)
    a2[1:2:end] = 3 * a1
    a2[2:2:end] = 5 * a1
    aold1 = rand(ndofs(dh1))
    aold2 = zeros(ndofs(dh2))
    aold2[1:2:end] = aold1;
    aold2[2:2:end] = aold1;

    for threading in (false, true)
        for autodiffbuffer in (false, true)
            for singledomain in (true, false)
                if singledomain    
                    d1 = setup_domainbuffer(DomainSpec(dh1, MA(), cvu); a = a1, threading, autodiffbuffer)
                    d2 = setup_domainbuffer(DomainSpec(dh2, MB(), cvv); a = a2, threading, autodiffbuffer)
                else
                    sets = Dict(k => getcellset(grid, k) for k in ("left", "right"))
                    d1 = setup_domainbuffers(Dict(k => DomainSpec(dh1, MA(), cvu; set) for (k, set) in sets); a = a1, threading, autodiffbuffer)
                    d2 = setup_domainbuffers(Dict(k => DomainSpec(dh2, MB(), cvv; set) for (k, set) in sets); a = a2, threading, autodiffbuffer)
                end
                # No setup-time coupling call: coupling is derived directly from whatever
                # `CoupledSimulations` is supplied to `work!`, fresh on every call.
                sim1 = Simulation(d1, a1, aold1)
                sim2 = Simulation(d2, a2, aold2)
                K = allocate_matrix(dh1)
                r = zeros(ndofs(dh1))
                assembler = start_assemble(K, r)
                Δt2 = 0.25
                set_time_increment!(d2, Δt2)
                work!(assembler, sim1, CoupledSimulations(b = sim2)) # Test

                # Changing the partner's time increment before the next staggered iteration
                # is picked up immediately: no persistent link to go stale (BUG-003 regression)
                Δt2 = 0.75
                set_time_increment!(d2, Δt2)
                assembler = start_assemble(K, r)
                work!(assembler, sim1, CoupledSimulations(b = sim2)) # Test

                # An independently replaced buffer (a genuinely different object from `d2`, as
                # occurs e.g. after `replace_material`) works transparently: there's no persistent
                # link that could go stale or mismatch, since coupling is derived fresh each call.
                d2_indep = FerriteAssembly.replace_material(d2, identity)
                sim2_indep = Simulation(d2_indep, a2, aold2)
                Δt2 = 0.4
                set_time_increment!(d2_indep, Δt2)
                assembler = start_assemble(K, r)
                work!(assembler, sim1, CoupledSimulations(b = sim2_indep)) # Test
            end
        end
    end

    # Task counts must match between primary and coupled partner: reusing a partner's per-task
    # buffers directly (no reallocation) is only race-free with a strict 1-to-1 mapping. A
    # sequential partner counts as having 1 task.
    d1_mismatch = setup_domainbuffer(DomainSpec(dh1, MA(), cvu); a = a1, threading = true, num_tasks = 3)
    for d2_mismatch in (
            setup_domainbuffer(DomainSpec(dh2, MB(), cvv); a = a2), # sequential = 1 task
            setup_domainbuffer(DomainSpec(dh2, MB(), cvv); a = a2, threading = true, num_tasks = 1),
        )
        sim1_mismatch = Simulation(d1_mismatch, a1, aold1)
        sim2_mismatch = Simulation(d2_mismatch, a2, aold2)
        Km = allocate_matrix(dh1)
        rm = zeros(ndofs(dh1))
        assemblerm = start_assemble(Km, rm)
        @test_throws ArgumentError work!(assemblerm, sim1_mismatch, CoupledSimulations(b = sim2_mismatch))
    end

    # Matching task counts (including a sequential partner matched to a 1-task primary) work,
    # reusing the partner's own buffer(s) directly - no reallocation.
    for (primary_num_tasks, partner_threading) in ((1, false), (2, true))
        d1_match = setup_domainbuffer(DomainSpec(dh1, MA(), cvu); a = a1, threading = true, num_tasks = primary_num_tasks)
        d2_match = partner_threading ?
            setup_domainbuffer(DomainSpec(dh2, MB(), cvv); a = a2, threading = true, num_tasks = primary_num_tasks) :
            setup_domainbuffer(DomainSpec(dh2, MB(), cvv); a = a2)
        sim1_match = Simulation(d1_match, a1, aold1)
        sim2_match = Simulation(d2_match, a2, aold2)
        Km = allocate_matrix(dh1)
        rm = zeros(ndofs(dh1))
        Δt2 = 0.6
        set_time_increment!(d2_match, Δt2)
        assemblerm = start_assemble(Km, rm)
        work!(assemblerm, sim1_match, CoupledSimulations(b = sim2_match)) # element_routine! asserts under real concurrency
    end
end

@testset "couple_buffers allocations" begin
    # Coupling should add O(1) allocation per `work!` call, not O(ncells): verify the extra
    # allocation from adding coupling doesn't scale with the number of cells.
    grid = generate_grid(Quadrilateral, (30, 30)) # 900 cells
    ip = Lagrange{RefQuadrilateral,1}()
    dh1 = close!(add!(DofHandler(grid), :u, ip))
    dh2 = close!(add!(DofHandler(grid), :v, ip))
    qr = QuadratureRule{RefQuadrilateral}(2)
    cv = CellValues(qr, ip, ip)
    struct MC end
    FerriteAssembly.element_routine!(Ke, re, state, ae, ::MC, cv, buffer) = fill!(Ke, 0)
    d1 = setup_domainbuffer(DomainSpec(dh1, MC(), cv))
    d2 = setup_domainbuffer(DomainSpec(dh2, MC(), cv))
    a1 = zeros(ndofs(dh1))
    a2 = zeros(ndofs(dh2))
    sim1 = Simulation(d1, a1, a1)
    sim2 = Simulation(d2, a2, a2)
    K = allocate_matrix(dh1)
    r = zeros(ndofs(dh1))

    assembler = start_assemble(K, r)
    work!(assembler, sim1) # compile/warmup, uncoupled
    assembler = start_assemble(K, r)
    work!(assembler, sim1, CoupledSimulations(b = sim2)) # compile/warmup, coupled

    assembler = start_assemble(K, r)
    nalloc_uncoupled = @allocated work!(assembler, sim1)
    assembler = start_assemble(K, r)
    nalloc_coupled = @allocated work!(assembler, sim1, CoupledSimulations(b = sim2))
    # 900 cells: any per-cell allocation (even tens of bytes) would show up as tens of KB here.
    @test (nalloc_coupled - nalloc_uncoupled) < 10_000
end

@testset "couple_buffers threaded allocations" begin
    # Threaded coupling reuses the partner's own per-task buffers directly (no reallocation), so
    # it should add ~no extra allocation over uncoupled work, and certainly not scale with ncells.
    ip = Lagrange{RefQuadrilateral,1}()
    qr = QuadratureRule{RefQuadrilateral}(2)
    cv = CellValues(qr, ip, ip)
    struct MD end
    FerriteAssembly.element_routine!(Ke, re, state, ae, ::MD, cv, buffer) = fill!(Ke, 0)

    function nalloc_work(ncells; coupled)
        # Coupling assumes matching grids, so both domains share the same grid/cell count.
        grid = generate_grid(Quadrilateral, (ncells, ncells))
        dh1 = close!(add!(DofHandler(grid), :u, ip))
        d1 = setup_domainbuffer(DomainSpec(dh1, MD(), cv); threading = true, num_tasks = 4)
        a1 = zeros(ndofs(dh1))
        sim1 = Simulation(d1, a1, a1)
        K = allocate_matrix(dh1)
        r = zeros(ndofs(dh1))
        coupled_sims = if coupled
            dh2 = close!(add!(DofHandler(grid), :v, ip))
            d2 = setup_domainbuffer(DomainSpec(dh2, MD(), cv); threading = true, num_tasks = 4)
            a2 = zeros(ndofs(dh2))
            CoupledSimulations(b = Simulation(d2, a2, a2))
        else
            CoupledSimulations()
        end
        assembler = start_assemble(K, r)
        work!(assembler, sim1, coupled_sims) # compile/warmup
        assembler = start_assemble(K, r)
        return @allocated work!(assembler, sim1, coupled_sims)
    end

    nalloc_uncoupled_30 = nalloc_work(30; coupled = false)
    nalloc_coupled_30 = nalloc_work(30; coupled = true)
    nalloc_coupled_60 = nalloc_work(60; coupled = true) # 4x more cells

    @test (nalloc_coupled_30 - nalloc_uncoupled_30) < 10_000 # reuse, not per-cell allocation
    @test nalloc_coupled_60 < 2 * nalloc_coupled_30 # doesn't scale with ncells

    # AutoDiffCellBuffer coupling: a per-task AutoDiffCellBuffer is (re)constructed each `work!`
    # call to relink `coupled_buffers` (a small, fixed-size cost - it does NOT recompute the
    # (expensive) JacobianConfig once the coupling structure is stable, and must not scale with
    # cell count). Reusing the *same* sim/buffers, unlike the plain-CellBuffer case above, so the
    # "type didn't change" fast path in `couple_buffers(::AutoDiffCellBuffer, ...)` is exercised.
    function nalloc_ad_coupled(ncells)
        grid = generate_grid(Quadrilateral, (ncells, ncells))
        dh1 = close!(add!(DofHandler(grid), :u, ip))
        dh2 = close!(add!(DofHandler(grid), :v, ip))
        d1 = setup_domainbuffer(DomainSpec(dh1, MD(), cv); threading = true, num_tasks = 4, autodiffbuffer = true)
        d2 = setup_domainbuffer(DomainSpec(dh2, MD(), cv); threading = true, num_tasks = 4, autodiffbuffer = true)
        a1 = zeros(ndofs(dh1))
        a2 = zeros(ndofs(dh2))
        sim1 = Simulation(d1, a1, a1)
        cs = CoupledSimulations(b = Simulation(d2, a2, a2))
        K = allocate_matrix(dh1)
        r = zeros(ndofs(dh1))
        assembler = start_assemble(K, r)
        work!(assembler, sim1, cs) # warmup call 1: coupled_buffers type changes empty -> coupled
        assembler = start_assemble(K, r)
        work!(assembler, sim1, cs) # warmup call 2: steady state (type already matches)
        assembler = start_assemble(K, r)
        n1 = @allocated work!(assembler, sim1, cs)
        assembler = start_assemble(K, r)
        n2 = @allocated work!(assembler, sim1, cs)
        return n1, n2
    end

    nalloc_ad_30a, nalloc_ad_30b = nalloc_ad_coupled(30)
    nalloc_ad_60a, _ = nalloc_ad_coupled(60) # 4x more cells

    @test nalloc_ad_30a == nalloc_ad_30b # steady state: no growth/repeated JacobianConfig rebuild
    @test nalloc_ad_60a < 2 * nalloc_ad_30a # doesn't scale with ncells
end
