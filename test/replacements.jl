@testset "replace_material" begin
    m_el = EE.LinearElastic(;E=1.0, ν=0.4)
    m_elx2 = EE.LinearElastic(;E=2.0, ν=0.4)
    f_repl1(::EE.LinearElastic) = m_elx2
    m_pl = EE.J2Plasticity(;E=1.0, ν=0.4, σ0=0.2, H=1.0)
    f_repl2(::EE.LinearElastic) = m_pl
    f_repl3(::EE.J2Plasticity) = m_el
    grid = generate_grid(Quadrilateral, (2,2))
    dh = DofHandler(grid); add!(dh, :u, 2); close!(dh)
    ip = Lagrange{2,RefCube,1}()
    qr = QuadratureRule{2,RefCube}(2)
    cv = CellVectorValues(qr, ip, ip)
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