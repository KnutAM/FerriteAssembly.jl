@testset "CoupledSimulations" begin
    # `@inferred` can't check constant propagation for property access (`.a` etc.): its macro
    # only accepts call expressions, and it infers based on the *runtime type* of arguments
    # passed to `getproperty`, not the literal property name baked into `x.a` syntax at the
    # call site (the case that actually matters, since that's how these are used everywhere).
    # `Base.return_types` on a closure containing the literal dot-access captures that.
    is_concrete_inferred(f, argtypes...) = begin
        rt = Base.return_types(f, argtypes)
        length(rt) == 1 && isconcretetype(rt[1])
    end

    grid = generate_grid(Quadrilateral, (2,2))
    addcellset!(grid, "left", x -> x[1] < eps())
    addcellset!(grid, "right", setdiff(1:getncells(grid), getcellset(grid, "left")))
    ip = Lagrange{RefQuadrilateral,1}()
    dh1 = close!(add!(DofHandler(grid), :u, ip))
    dh2 = close!(add!(DofHandler(grid), :v, ip^2))
    qr = QuadratureRule{RefQuadrilateral}(2)
    cvu = CellValues(qr, ip, ip)
    cvv = CellValues(qr, ip^2, ip)

    struct CS_MA end
    struct CS_MB end
    struct CS_MB2 end # distinct material type, used to exercise replace_material changing type
    # aold will be same in both cases (for both components in the case of MB)
    # a will be 3 times larger for first component in MB, and 5 times for second component
    # State will be 6 times larger for MB, obtained by multiplying the function values by factor 2
    FerriteAssembly.create_cell_state(::CS_MA, cv, x, ae, args...) = [function_value(cv, i, ae) for i in 1:getnquadpoints(cv)]
    FerriteAssembly.create_cell_state(::CS_MB, cv, x, ae, args...) = [2 * function_value(cv, i, ae)[1] for i in 1:getnquadpoints(cv)]
    FerriteAssembly.create_cell_state(::CS_MB2, cv, x, ae, args...) = nothing

    # Set to a Float64 (not NaN) by the BUG-003 regression below to check that the partner's
    # (`:b`'s) *own* current time increment is observed, independent of the reader's own Δt.
    expected_b_dt = Ref(NaN)
    expected_b_material = Ref{DataType}(CS_MB)
    function FerriteAssembly.element_routine!(Ke, re, state, ae, m::CS_MA, cv, buffer)
        cb_b = FerriteAssembly.get_coupled_buffer(buffer, :b)
        @test FerriteAssembly.get_material(cb_b) isa expected_b_material[]
        if expected_b_material[] === CS_MB
            @test 3 * ae ≈ FerriteAssembly.get_ae(cb_b)[1:2:end]
            @test 5 * ae ≈ FerriteAssembly.get_ae(cb_b)[2:2:end]
            @test FerriteAssembly.get_aeold(buffer) ≈ FerriteAssembly.get_aeold(cb_b)[1:2:end]
            @test FerriteAssembly.get_aeold(buffer) ≈ FerriteAssembly.get_aeold(cb_b)[2:2:end]
            @test 6 * state ≈ FerriteAssembly.get_state(cb_b)
        end
        isnan(expected_b_dt[]) || @test FerriteAssembly.get_time_increment(cb_b) == expected_b_dt[]
        # Present only in the 3-member mutual-coupling test below; checks that :b and :c are
        # not positionally swapped when a reader has two distinct partners.
        if haskey(FerriteAssembly.get_coupled_buffers(buffer), :c)
            cb_c = FerriteAssembly.get_coupled_buffer(buffer, :c)
            @test 7 * ae ≈ FerriteAssembly.get_ae(cb_c)
        end
    end
    function FerriteAssembly.element_routine!(Ke, re, state, ae, m::CS_MB, cv, buffer)
        nothing # Only assembled from `:a`'s perspective in these tests
    end
    function FerriteAssembly.element_routine!(Ke, re, state, ae, m::CS_MB2, cv, buffer)
        nothing
    end

    a1 = rand(ndofs(dh1))
    a2 = zeros(ndofs(dh2))
    @assert length(a1) * 2 == length(a2)
    a2[1:2:end] = 3 * a1
    a2[2:2:end] = 5 * a1
    aold1 = rand(ndofs(dh1))
    aold2 = zeros(ndofs(dh2))
    aold2[1:2:end] = aold1
    aold2[2:2:end] = aold1

    @testset "threading=$threading, autodiffbuffer=$autodiffbuffer, singledomain=$singledomain" for
            threading in (false, true), autodiffbuffer in (false, true), singledomain in (true, false)
        if singledomain
            d1 = setup_domainbuffer(DomainSpec(dh1, CS_MA(), cvu); a = a1, threading, autodiffbuffer)
            d2 = setup_domainbuffer(DomainSpec(dh2, CS_MB(), cvv); a = a2, threading, autodiffbuffer)
        else
            sets = Dict(k => getcellset(grid, k) for k in ("left", "right"))
            d1 = setup_domainbuffers(Dict(k => DomainSpec(dh1, CS_MA(), cvu; set) for (k, set) in sets); a = a1, threading, autodiffbuffer)
            d2 = setup_domainbuffers(Dict(k => DomainSpec(dh2, CS_MB(), cvv; set) for (k, set) in sets); a = a2, threading, autodiffbuffer)
        end
        sim1 = Simulation(d1, a1, aold1)
        sim2 = Simulation(d2, a2, aold2)
        g = CoupledSimulations((a = sim1,); refs = (b = sim2,))
        @test g.a isa FerriteAssembly.CoupledSimulation
        @test g.b === sim2 # refs are the plain source Simulation

        # A CoupledSimulation member forwards ordinary Simulation property access: `.a`/
        # `.aold` are the same global vectors (shared by reference, not copied), `.db` is the
        # *rebuilt* (coupled) domain buffer, not the original source `d1`.
        @test g.a.a === a1
        @test g.a.aold === aold1
        @test g.a.db === g.a.sim.db
        @test g.a.db !== d1

        # Property access must constant-propagate to a single concrete type (not a Union
        # across every forwarding branch) for both the member handle's own getproperty
        # override and the group's.
        @test is_concrete_inferred(csim -> csim.a, typeof(g.a))
        @test is_concrete_inferred(csim -> csim.aold, typeof(g.a))
        @test is_concrete_inferred(csim -> csim.db, typeof(g.a))
        @test is_concrete_inferred(csim -> csim.sim, typeof(g.a))
        @test is_concrete_inferred(csim -> csim.partners, typeof(g.a))
        @test is_concrete_inferred(grp -> grp.a, typeof(g))

        # Forwarded properties must tab-complete, not just the two real struct fields.
        propnames = propertynames(g.a)
        @test :a in propnames
        @test :aold in propnames
        @test :db in propnames
        @test :sim in propnames
        @test :partners in propnames

        K = allocate_matrix(dh1)
        r = zeros(ndofs(dh1))
        assembler = start_assemble(K, r)
        expected_b_dt[] = NaN
        work!(assembler, g.a) # Runs the @test's inside element_routine!

        # BUG-003 regression: changing a ref's Δt between two work! calls must be observed
        # by every reader task, not just task 1, without re-working the ref itself.
        set_time_increment!(g.b, 1.23)
        expected_b_dt[] = 1.23
        work!(assembler, g.a)
        set_time_increment!(g.b, 4.56)
        expected_b_dt[] = 4.56
        work!(assembler, g.a) # element_routine! asserts Δt equality itself
        expected_b_dt[] = NaN

        # Stable buffer/config identity across repeated work! calls (no rebuild per call)
        get_ib() = singledomain ? FerriteAssembly.get_itembuffer(g.a) : FerriteAssembly.get_itembuffer(g.a, "left")
        ib1 = FerriteAssembly.get_base(get_ib())
        work!(assembler, g.a)
        ib2 = FerriteAssembly.get_base(get_ib())
        @test ib1 === ib2

        # Coupling does not scale allocations with cell count (allow generous fixed overhead
        # for task-spawn/chunk machinery; this is a smoke check, not a scaling sweep).
        work!(assembler, g.a) # warm up (compile)
        work!(assembler, g.a) # warm up again to be safe against any first-use effects
        nalloc = @allocated work!(assembler, g.a)
        @test nalloc < 2_000_000
    end

    @testset "coupling allocations do not scale with cell count" begin
        # The single-mesh check above only bounds allocations against a fixed ceiling on one
        # 4-cell grid; it cannot detect a small per-cell allocation (e.g. a reintroduced
        # per-cell wrapper/config construction) that would still be far below that ceiling.
        # Compare a much larger mesh against a tiny one instead: coupling-specific overhead
        # (wrapper/config construction, task-spawn/chunk machinery) is paid once per `work!`
        # call, not per cell, so it must not grow materially with cell count.
        struct CS_AllocA end
        struct CS_AllocB end
        struct CS_AllocA0 end # uncoupled baseline: same per-cell work, no partner access
        FerriteAssembly.create_cell_state(::CS_AllocA, args...) = nothing
        FerriteAssembly.create_cell_state(::CS_AllocB, args...) = nothing
        FerriteAssembly.create_cell_state(::CS_AllocA0, args...) = nothing
        function FerriteAssembly.element_routine!(Ke, re, state, ae, ::CS_AllocA, cv, buffer)
            cb = FerriteAssembly.get_coupled_buffer(buffer, :b)
            ae_p = FerriteAssembly.get_ae(cb)
            @inbounds for i in eachindex(re)
                re[i] += ae_p[i]
            end
            return nothing
        end
        FerriteAssembly.element_routine!(Ke, re, state, ae, ::CS_AllocB, cv, buffer) = nothing
        function FerriteAssembly.element_routine!(Ke, re, state, ae, ::CS_AllocA0, cv, buffer)
            @inbounds for i in eachindex(re)
                re[i] += ae[i]
            end
            return nothing
        end

        function build_alloc_group(n)
            grid_ = generate_grid(Quadrilateral, (n, n))
            ip_ = Lagrange{RefQuadrilateral,1}()
            dhA = close!(add!(DofHandler(grid_), :u, ip_))
            dhB = close!(add!(DofHandler(grid_), :v, ip_^2))
            cvA = CellValues(qr, ip_, ip_)
            cvB = CellValues(qr, ip_^2, ip_)
            aA = zeros(ndofs(dhA))
            aB = zeros(ndofs(dhB))
            dA = setup_domainbuffer(DomainSpec(dhA, CS_AllocA(), cvA); a = aA)
            dB = setup_domainbuffer(DomainSpec(dhB, CS_AllocB(), cvB); a = aB)
            simA = Simulation(dA, aA, zeros(ndofs(dhA)))
            simB = Simulation(dB, aB, zeros(ndofs(dhB)))
            return CoupledSimulations((a = simA,); refs = (b = simB,)), dhA
        end
        function build_baseline(n)
            grid_ = generate_grid(Quadrilateral, (n, n))
            ip_ = Lagrange{RefQuadrilateral,1}()
            dhA = close!(add!(DofHandler(grid_), :u, ip_))
            cvA = CellValues(qr, ip_, ip_)
            aA = zeros(ndofs(dhA))
            dA = setup_domainbuffer(DomainSpec(dhA, CS_AllocA0(), cvA); a = aA)
            return Simulation(dA, aA, zeros(ndofs(dhA))), dhA
        end
        function measure_alloc(sim_or_group, dh_)
            K = allocate_matrix(dh_)
            r = zeros(ndofs(dh_))
            asm = start_assemble(K, r)
            work!(asm, sim_or_group) # warm up (compile)
            work!(asm, sim_or_group) # warm up again
            return @allocated work!(asm, sim_or_group)
        end
        g_small, dh_small = build_alloc_group(2)
        g_large, dh_large = build_alloc_group(20) # 100x the cells of g_small
        nalloc_small = measure_alloc(g_small.a, dh_small)
        nalloc_large = measure_alloc(g_large.a, dh_large)

        # Coupling-specific overhead relative to an uncoupled baseline doing equivalent
        # per-cell work: both must be exactly zero (ordinary sequential assembly is
        # allocation-free), so even a small per-cell allocation reintroduced by coupling
        # would be caught, not just growth that outpaces cell count.
        base_small = measure_alloc(build_baseline(2)...)
        base_large = measure_alloc(build_baseline(20)...)
        @test nalloc_small - base_small == 0
        @test nalloc_large - base_large == 0
    end

    @testset "threaded reader (1 task) with sequential partner" begin
        # validate_domain_pair explicitly allows this (a sequential partner counts as 1 slot,
        # matching a threaded reader with exactly 1 task); work! must not throw when scattering
        # partners before dispatch, even though the partner has no task-local buffers to
        # scatter into.
        expected_b_dt[] = NaN
        expected_b_material[] = CS_MB
        d1 = setup_domainbuffer(DomainSpec(dh1, CS_MA(), cvu); a = a1, threading = true, num_tasks = 1)
        d2 = setup_domainbuffer(DomainSpec(dh2, CS_MB(), cvv); a = a2, threading = false)
        sim1 = Simulation(d1, a1, aold1)
        sim2 = Simulation(d2, a2, aold2)
        g = CoupledSimulations((a = sim1,); refs = (b = sim2,))
        K = allocate_matrix(dh1)
        r = zeros(ndofs(dh1))
        work!(start_assemble(K, r), g.a)
    end

    @testset "mutual coupling (3 members)" begin
        struct CS_MC end
        FerriteAssembly.create_cell_state(::CS_MC, cv, x, ae, args...) = nothing
        function FerriteAssembly.element_routine!(Ke, re, state, ae, ::CS_MC, cv, buffer)
            nothing
        end
        ip3 = Lagrange{RefQuadrilateral,1}()
        dh3 = close!(add!(DofHandler(grid), :w, ip3))
        cv3 = CellValues(qr, ip3, ip3)
        # Same dof ordering as dh1 (same grid/interpolation/single scalar field), so a
        # component-wise multiple of a1 gives an independently checkable per-cell value,
        # matching the existing a2/a1 pattern used for the :b partner above.
        a3 = 7 * a1
        aold3 = zeros(ndofs(dh3))
        expected_b_material[] = CS_MB
        for threading in (false, true)
            d1 = setup_domainbuffer(DomainSpec(dh1, CS_MA(), cvu); a = a1, threading)
            d2 = setup_domainbuffer(DomainSpec(dh2, CS_MB(), cvv); a = a2, threading)
            d3 = setup_domainbuffer(DomainSpec(dh3, CS_MC(), cv3); a = a3, threading)
            sim1 = Simulation(d1, a1, aold1)
            sim2 = Simulation(d2, a2, aold2)
            sim3 = Simulation(d3, a3, aold3)
            g = CoupledSimulations((a = sim1, b = sim2, c = sim3))
            @test g.a isa FerriteAssembly.CoupledSimulation
            @test g.b isa FerriteAssembly.CoupledSimulation
            @test g.c isa FerriteAssembly.CoupledSimulation
            cb_b = FerriteAssembly.get_coupled_buffers(FerriteAssembly.get_base(FerriteAssembly.get_itembuffer(g.a)))
            @test haskey(cb_b, :b) && haskey(cb_b, :c)
            # b and c views from a are nonrecursive: plain CellBuffer, no further coupling
            @test !hasmethod(FerriteAssembly.get_coupled_buffers, Tuple{typeof(cb_b.b)})

            # Actually work! the reader with two distinct partners: CS_MA's element_routine!
            # checks both :b (2x/2y-scaled dof/state values) and :c (7x-scaled dof values),
            # so a positional mix-up between the two partner NamedTuples (buffers vs.
            # simulations) would fail here even though it could pass a construction-only check.
            expected_b_dt[] = NaN
            K = allocate_matrix(dh1)
            r = zeros(ndofs(dh1))
            work!(start_assemble(K, r), g.a)
        end
    end

    @testset "autodiff through coupling: numerical agreement" begin
        # `CS_AD_Reader` defines only `element_residual!` (no `element_routine!`), so its Ke
        # is genuinely computed via ForwardDiff through the coupled partner buffer, not a
        # hand-written Ke. `CS_AD_Reader_manual` computes the analytically-known Ke (c*I,
        # since the residual is elementwise linear) directly via `element_routine!`, coupled
        # to the same partner. Assembling both through the same mesh/assembly machinery and
        # comparing the results verifies the AD-through-coupling path numerically, for both
        # sequential and threaded execution.
        struct CS_AD_Reader
            c::Float64
        end
        struct CS_AD_Reader_manual
            c::Float64
        end
        struct CS_AD_Partner end
        FerriteAssembly.create_cell_state(::CS_AD_Reader, args...) = nothing
        FerriteAssembly.create_cell_state(::CS_AD_Reader_manual, args...) = nothing
        FerriteAssembly.create_cell_state(::CS_AD_Partner, args...) = nothing

        function FerriteAssembly.element_residual!(re, state, ae, m::CS_AD_Reader, cv, buffer)
            ae_p = FerriteAssembly.get_ae(FerriteAssembly.get_coupled_buffer(buffer, :p))
            re .= m.c .* ae .- ae_p
            return nothing
        end
        function FerriteAssembly.element_routine!(Ke, re, state, ae, m::CS_AD_Reader_manual, cv, buffer)
            ae_p = FerriteAssembly.get_ae(FerriteAssembly.get_coupled_buffer(buffer, :p))
            re .= m.c .* ae .- ae_p
            fill!(Ke, 0)
            for i in axes(Ke, 1)
                Ke[i, i] = m.c
            end
            return nothing
        end
        function FerriteAssembly.element_routine!(Ke, re, state, ae, ::CS_AD_Partner, cv, buffer)
            nothing
        end

        ipr = Lagrange{RefQuadrilateral,1}()
        dhr = close!(add!(DofHandler(grid), :r, ipr))
        cvr = CellValues(qr, ipr, ipr)
        ar = rand(ndofs(dhr))
        ap = rand(ndofs(dhr))
        aold_dummy = zeros(ndofs(dhr))
        c = 2.5

        for threading in (false, true)
            dr_ad = setup_domainbuffer(DomainSpec(dhr, CS_AD_Reader(c), cvr); a = ar, threading, autodiffbuffer=true)
            dr_man = setup_domainbuffer(DomainSpec(dhr, CS_AD_Reader_manual(c), cvr); a = ar, threading)
            dp = setup_domainbuffer(DomainSpec(dhr, CS_AD_Partner(), cvr); a = ap, threading)
            simr_ad = Simulation(dr_ad, ar, aold_dummy)
            simr_man = Simulation(dr_man, ar, aold_dummy)
            simp = Simulation(dp, ap, aold_dummy)
            g_ad = CoupledSimulations((r = simr_ad,); refs = (p = simp,))
            g_man = CoupledSimulations((r = simr_man,); refs = (p = simp,))
            K_ad = allocate_matrix(dhr); r_ad = zeros(ndofs(dhr))
            K_man = allocate_matrix(dhr); r_man = zeros(ndofs(dhr))
            work!(start_assemble(K_ad, r_ad), g_ad.r)
            work!(start_assemble(K_man, r_man), g_man.r)
            @test Matrix(K_ad) ≈ Matrix(K_man)
            @test r_ad ≈ r_man
        end
    end

    @testset "replace_material through group" begin
        d1 = setup_domainbuffer(DomainSpec(dh1, CS_MA(), cvu); a = a1)
        d2 = setup_domainbuffer(DomainSpec(dh2, CS_MB(), cvv); a = a2)
        sim1 = Simulation(d1, a1, aold1)
        sim2 = Simulation(d2, a2, aold2)
        g = CoupledSimulations((a = sim1,); refs = (b = sim2,))
        K = allocate_matrix(dh1); r = zeros(ndofs(dh1))
        expected_b_dt[] = NaN
        expected_b_material[] = CS_MB
        work!(start_assemble(K, r), g.a) # exercise before replacement, observes CS_MB

        g2 = FerriteAssembly.replace_material(g, :b, m -> CS_MB2()) # changes the material TYPE
        @test g2.a isa FerriteAssembly.CoupledSimulation
        @test FerriteAssembly.get_material(g2.b) isa CS_MB2
        @test FerriteAssembly.get_material(g.b) isa CS_MB # old group/handle untouched
        @test g2.a !== g.a

        # New handle actually works and observes the new material through coupling
        expected_b_material[] = CS_MB2
        work!(start_assemble(K, r), g2.a)
        # Old handle, worked again, still observes the original material (not silently switched)
        expected_b_material[] = CS_MB
        work!(start_assemble(K, r), g.a)
        expected_b_material[] = CS_MB

        @test_throws ArgumentError FerriteAssembly.replace_material(g, :nope, identity)
        @test_throws ArgumentError FerriteAssembly.replace_material(g.a, identity)
    end

    @testset "validation errors" begin
        d1 = setup_domainbuffer(DomainSpec(dh1, CS_MA(), cvu); a = a1)
        d2 = setup_domainbuffer(DomainSpec(dh2, CS_MB(), cvv); a = a2)
        sim1 = Simulation(d1, a1, aold1)
        sim2 = Simulation(d2, a2, aold2)

        @test_throws ArgumentError CoupledSimulations(NamedTuple()) # empty primaries
        @test_throws ArgumentError CoupledSimulations((a = sim1,); refs = (a = sim2,)) # duplicate name

        # Different grid
        grid2 = generate_grid(Quadrilateral, (2,2))
        ip2 = Lagrange{RefQuadrilateral,1}()
        dh2b = close!(add!(DofHandler(grid2), :v, ip2^2))
        a2b = zeros(ndofs(dh2b))
        d2b = setup_domainbuffer(DomainSpec(dh2b, CS_MB(), CellValues(qr, ip2^2, ip2)); a = a2b)
        sim2b = Simulation(d2b, a2b, zeros(ndofs(dh2b)))
        @test_throws ArgumentError CoupledSimulations((a = sim1,); refs = (b = sim2b,))

        # Missing domain coverage (multi-domain reader, partner missing a domain)
        sets = Dict(k => getcellset(grid, k) for k in ("left", "right"))
        d1m = setup_domainbuffers(Dict(k => DomainSpec(dh1, CS_MA(), cvu; set) for (k, set) in sets); a = a1)
        d2m_partial = Dict("left" => setup_domainbuffer(DomainSpec(dh2, CS_MB(), cvv; set=sets["left"]); a = a2))
        sim1m = Simulation(d1m, a1, aold1)
        sim2m_partial = Simulation(d2m_partial, a2, aold2)
        @test_throws ArgumentError CoupledSimulations((a = sim1m,); refs = (b = sim2m_partial,))

        # Mixed single/dictionary coupling
        @test_throws ArgumentError CoupledSimulations((a = sim1m,); refs = (b = sim2,))
        @test_throws ArgumentError CoupledSimulations((a = sim1,); refs = (b = Simulation(d2m_partial, a2, aold2),))

        # Incompatible task counts
        d1t = setup_domainbuffer(DomainSpec(dh1, CS_MA(), cvu); a = a1, threading=true, num_tasks=2)
        d2t = setup_domainbuffer(DomainSpec(dh2, CS_MB(), cvv); a = a2, threading=true, num_tasks=3)
        sim1t = Simulation(d1t, a1, aold1)
        sim2t = Simulation(d2t, a2, aold2)
        @test_throws ArgumentError CoupledSimulations((a = sim1t,); refs = (b = sim2t,))

        # Duplicate storage: same domain buffer object used under two member names
        @test_throws ArgumentError CoupledSimulations((a = sim1,); refs = (b = sim2, c = sim2))

        # Unsupported buffer kind (facet buffers): rejected with an actionable ArgumentError,
        # not a MethodError, both as a partner and as a sole primary with no partners at all
        # (storage-identity validation runs for every member, regardless of pairing).
        struct CS_MFacet end
        dh_f = close!(add!(DofHandler(grid), :f, ip))
        fv = FacetValues(FacetQuadratureRule{RefQuadrilateral}(2), ip)
        d_facet = setup_domainbuffer(DomainSpec(dh_f, CS_MFacet(), fv; set=getfacetset(grid, "left")))
        a_f = zeros(ndofs(dh_f))
        sim_facet = Simulation(d_facet, a_f, zeros(ndofs(dh_f)))
        @test_throws ArgumentError CoupledSimulations((a = sim_facet,))
        @test_throws ArgumentError CoupledSimulations((a = sim1,); refs = (b = sim_facet,))
    end
end
