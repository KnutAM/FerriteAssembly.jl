@testset "Item buffers" begin
    struct IB_Material end
    struct IB_Partner end
    struct IB_FacetMaterial end
    FerriteAssembly.allocate_cell_cache(::IB_Material, ::Any) = [1.0]
    FerriteAssembly.allocate_facet_cache(::IB_FacetMaterial, ::Any) = [2.0]

    grid = generate_grid(Quadrilateral, (2, 2))
    ip = Lagrange{RefQuadrilateral, 1}()
    dh = close!(add!(add!(DofHandler(grid), :u, ip^2), :p, ip))
    sdh = dh.subdofhandlers[1]
    qr = QuadratureRule{RefQuadrilateral}(2)
    cv = (u = CellValues(qr, ip^2, ip), p = CellValues(qr, ip, ip))
    fv = (u = FacetValues(FacetQuadratureRule{RefQuadrilateral}(2), ip^2, ip),
          p = FacetValues(FacetQuadratureRule{RefQuadrilateral}(2), ip, ip))
    a = rand(ndofs(dh))
    aold = rand(ndofs(dh))
    cellnum = 3

    @testset "DofBuffer" begin
        d = FA.DofBuffer(sdh)
        @test d.dofrange === (u = dof_range(sdh, :u), p = dof_range(sdh, :p))
        @test all(length(getfield(d, name)) == ndofs_per_cell(sdh) for name in (:ae_old, :ae, :re, :dofs))
        @test size(d.Ke) == (ndofs_per_cell(sdh), ndofs_per_cell(sdh))

        FA.set_dofs!(d, dh, cellnum)
        @test d.dofs == celldofs(dh, cellnum)

        db = setup_domainbuffer(DomainSpec(dh, nothing, cv))
        fill!(d.Ke, 1); fill!(d.re, 1); fill!(d.dofs, 0)
        FA.reinit_buffer!(d, Simulation(db, a, nothing), cellnum)
        @test d.dofs == celldofs(dh, cellnum)
        @test d.ae == a[celldofs(dh, cellnum)]
        @test all(isnan, d.ae_old)
        @test all(iszero, d.Ke) && all(iszero, d.re)
        FA.reinit_buffer!(d, Simulation(db, nothing, aold), cellnum)
        @test all(isnan, d.ae)
        @test d.ae_old == aold[celldofs(dh, cellnum)]

        d_local = FA.create_local(d)
        for name in (:ae_old, :ae, :re, :Ke, :dofs)
            @test isequal(getfield(d_local, name), getfield(d, name))
            @test getfield(d_local, name) !== getfield(d, name)
        end
        @test d_local.dofrange === d.dofrange
    end

    # The content of the buffer for the item belonging to cell number `cellnum`
    function test_item_content(buffer, dofbuffer)
        dofs = celldofs(dh, cellnum)
        @test FA.get_dofbuffer(buffer) === dofbuffer
        @test celldofs(buffer) === dofbuffer.dofs
        @test celldofs(buffer) == dofs
        @test FA.get_ae(buffer) === dofbuffer.ae
        @test FA.get_ae(buffer) == a[dofs]
        @test FA.get_aeold(buffer) === dofbuffer.ae_old
        @test FA.get_aeold(buffer) == aold[dofs]
        @test FA.get_re(buffer) === dofbuffer.re
        @test FA.get_Ke(buffer) === dofbuffer.Ke
        @test all(iszero, FA.get_re(buffer)) && all(iszero, FA.get_Ke(buffer))
        @test Ferrite.getfieldnames(buffer) === (:u, :p)
        get_range_u(b) = dof_range(b, :u)
        @test (@inferred get_range_u(buffer)) === dof_range(sdh, :u)
        @test dof_range(buffer, :p) === dof_range(sdh, :p)
        @test getcoordinates(buffer) == getcoordinates(grid, cellnum)
        @test cellid(buffer) == cellnum
    end

    @testset "CellBuffer" begin
        db = setup_domainbuffer(DomainSpec(dh, IB_Material(), cv; user_data = :cell_data))
        cb = FA.get_itembuffer(db)
        @test cb isa FA.CellBuffer
        fill!(FA.get_Ke(cb), 1); fill!(FA.get_re(cb), 1)
        FA.reinit_buffer!(cb, Simulation(db, a, aold), cellnum)
        test_item_content(cb, cb.dofbuffer)
        @test shape_value(FA.get_values(cb).u, 1, 1) == shape_value(cv.u, 1, 1)
        @test FA.get_user_data(cb) === :cell_data
        @test FA.get_user_cache(cb) == [1.0]

        cb_local = FA.create_local(cb)
        @test cb_local.dofbuffer !== cb.dofbuffer
        @test FA.get_ae(cb_local) == FA.get_ae(cb)
        @test FA.get_ae(cb_local) !== FA.get_ae(cb)
        @test celldofs(cb_local) !== celldofs(cb)
        @test FA.get_user_data(cb_local) === FA.get_user_data(cb)
        @test FA.get_user_cache(cb_local) !== FA.get_user_cache(cb)
    end

    @testset "FacetBuffer" begin
        set = [FacetIndex(cellnum, 4), FacetIndex(1, 4)]
        db = setup_domainbuffer(DomainSpec(dh, IB_FacetMaterial(), fv; set, user_data = :facet_data))
        fb = FA.get_itembuffer(db)
        @test fb isa FA.FacetBuffer
        fill!(FA.get_Ke(fb), 1); fill!(FA.get_re(fb), 1)
        FA.reinit_buffer!(fb, Simulation(db, a, aold), FacetIndex(cellnum, 4))
        test_item_content(fb, fb.dofbuffer)
        @test Ferrite.getcurrentfacet(FA.get_values(fb).u) == 4
        @test FA.get_user_data(fb) === :facet_data
        @test FA.get_user_cache(fb) == [2.0]

        fb_local = FA.create_local(fb)
        @test fb_local.dofbuffer !== fb.dofbuffer
        @test FA.get_ae(fb_local) == FA.get_ae(fb)
        @test FA.get_ae(fb_local) !== FA.get_ae(fb)
        @test FA.get_user_cache(fb_local) !== FA.get_user_cache(fb)
    end

    @testset "Wrapped cell buffers" begin
        setup_sim(material; kwargs...) = Simulation(setup_domainbuffer(DomainSpec(dh, material, cv; user_data = :cell_data); kwargs...), a, aold)
        plain_sim = setup_sim(IB_Material())
        ad_sim = setup_sim(IB_Material(); autodiffbuffer = true)
        coupled = CoupledSimulations((main = setup_sim(IB_Material()), partner = setup_sim(IB_Partner())))
        ad_coupled = CoupledSimulations((main = setup_sim(IB_Material(); autodiffbuffer = true), partner = setup_sim(IB_Partner())))

        plain_buffer = FA.get_itembuffer(plain_sim)
        ad_buffer = FA.get_itembuffer(ad_sim)
        coupled_buffer = FA.get_itembuffer(coupled.main)
        ad_coupled_buffer = FA.get_itembuffer(ad_coupled.main)
        @test ad_buffer isa FA.AutoDiffCellBuffer{<:FA.CellBuffer}
        @test coupled_buffer isa FA.CoupledCellBuffer
        @test ad_coupled_buffer isa FA.AutoDiffCellBuffer{<:FA.CoupledCellBuffer}

        identical_accessors = (FA.get_dofbuffer, FA.get_ae, FA.get_aeold, FA.get_re, FA.get_Ke,
            FA.get_material, FA.get_values, FA.get_state, FA.get_old_state, FA.get_user_data,
            FA.get_user_cache, celldofs, getcoordinates, Ferrite.getfieldnames)
        @testset "$(nameof(typeof(buffer))){$(nameof(typeof(FA.get_wrapped(buffer))))}" for (buffer, sim) in (
                (ad_buffer, ad_sim), (coupled_buffer, coupled.main), (ad_coupled_buffer, ad_coupled.main))
            @test buffer isa FA.AbstractCellBufferWrapper
            wrapped = FA.get_wrapped(buffer)
            FA.reinit_buffer!(buffer, sim, cellnum)
            test_item_content(buffer, FA.get_dofbuffer(wrapped))
            for accessor in identical_accessors
                @test accessor(buffer) === accessor(wrapped)
            end
            @test cellid(buffer) === cellid(wrapped)
            @test dof_range(buffer, :p) === dof_range(wrapped, :p)
            @test FA.get_user_data(buffer) === :cell_data
            @test FA.get_user_cache(buffer) == [1.0]

            FA.set_time_increment!(buffer, 0.25)
            @test FA.get_time_increment(buffer) == 0.25
            @test FA.get_time_increment(wrapped) == 0.25
        end

        # Only coupled buffers have coupled buffers
        @test !hasmethod(FA.get_coupled_buffers, Tuple{typeof(plain_buffer)})
        @test_throws MethodError FA.get_coupled_buffers(ad_buffer)
        partner_buffers = FA.get_coupled_buffers(coupled_buffer)
        @test partner_buffers === coupled_buffer.partner_buffers
        @test keys(partner_buffers) == (:partner,)
        @test FA.get_coupled_buffer(coupled_buffer, :partner) === partner_buffers.partner
        @test FA.get_coupled_buffers(ad_coupled_buffer) === FA.get_coupled_buffers(FA.get_wrapped(ad_coupled_buffer))
        @test FA.get_material(FA.get_coupled_buffer(ad_coupled_buffer, :partner)) isa IB_Partner

        # A coupled buffer requires a coupled simulation
        @test_throws MethodError FA.reinit_buffer!(coupled_buffer, plain_sim, cellnum)
    end
end
