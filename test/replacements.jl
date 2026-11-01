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

    # Domain-selective replacement
    bs3 = FerriteAssembly.replace_material(buffers, "a", f_repl1)
    @test FerriteAssembly.get_material(bs3, "a") === m_elx2
    @test FerriteAssembly.get_material(bs3, "b") === m_pl # unchanged, f_repl1 not applied here
    @test bs3["b"] === buffers["b"] # copied by reference
    @test_throws ArgumentError FerriteAssembly.replace_material(buffers, "c", f_repl1)

    @testset "replace_material shares item-buffer scratch (not just material) by reference" begin
        # This is the exact mechanism that makes `replace_material`-derived siblings - or a
        # sibling and its original template - unsafe to `work!` concurrently: everything
        # except `material` itself (ae, re, Ke, dofs, coords, cellvalues, user_data,
        # user_cache, and the buffer's states) is the *same* object, shared by reference.
        # See the "Concurrent work! calls" section in docs/src/DomainBuffers/Setup.md.
        b0 = setup_domainbuffer(dspec; threading=true, num_tasks=2)
        b1 = FerriteAssembly.replace_material(b0, f_repl1)

        @test FerriteAssembly.get_material(b0) !== FerriteAssembly.get_material(b1) # the one thing that changes
        @test FerriteAssembly.get_state(b0) === FerriteAssembly.get_state(b1)
        @test FerriteAssembly.get_old_state(b0) === FerriteAssembly.get_old_state(b1)

        ib0, ib1 = FerriteAssembly.get_itembuffer(b0), FerriteAssembly.get_itembuffer(b1)
        @test FerriteAssembly.get_base(ib0) !== FerriteAssembly.get_base(ib1) # new wrapper object...
        @test FerriteAssembly.get_base(ib0).ae === FerriteAssembly.get_base(ib1).ae # ...around the same scratch
        for i in 1:2
            cb0, cb1 = FerriteAssembly.get_local(ib0, i), FerriteAssembly.get_local(ib1, i)
            @test cb0 !== cb1
            @test cb0.material !== cb1.material
            for fname in (:ae_old, :ae, :re, :Ke, :dofs, :coords, :cellvalues, :user_data, :user_cache)
                @test getfield(cb0, fname) === getfield(cb1, fname)
            end
        end
    end
end
