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
    grid = generate_grid(Quadrilateral, (1,2))
    ip = Lagrange{RefQuadrilateral,1}()
    dh1 = close!(add!(DofHandler(grid), :u, ip))
    dh2 = close!(add!(DofHandler(grid), :v, ip^2))
    qr = QuadratureRule{RefQuadrilateral}(1)
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

    # Test case to check that values have been updated correctly
    function FerriteAssembly.element_routine!(Ke, re, state, ae, m::MA, cv, buffer)
        cb_b = buffer.coupled_buffers.b
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
    end
    
    # Simple domains, not threaded
    a1 = rand(ndofs(dh1))
    a2 = zeros(ndofs(dh2))
    @assert length(a1) * 2 == length(a2)
    a2[1:2:end] = 3 * a1
    a2[2:2:end] = 5 * a1
    aold1 = rand(ndofs(dh1))
    aold2 = zeros(ndofs(dh2))
    aold2[1:2:end] = aold1;
    aold2[2:2:end] = aold1;

    d1 = setup_domainbuffer(DomainSpec(dh1, MA(), cvu); a = a1)
    d2 = setup_domainbuffer(DomainSpec(dh2, MB(), cvv); a = a2)
    d1 = couple_buffers(d1; b = d2)
    sim1 = Simulation(d1, a1, aold1)
    sim2 = Simulation(d2, a2, aold2)
    K = allocate_matrix(dh1)
    r = zeros(ndofs(dh1))
    assembler = start_assemble(K, r)
    work!(assembler, sim1, CoupledSimulations(b = sim2)) # Test
end
