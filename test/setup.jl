@testset "setup" begin
    
    @testset "cells" begin
        struct CellMatCache
            dummy::Vector{Int}
        end
        CellMatCache() = CellMatCache([1])
        FerriteAssembly.allocate_cell_cache(::CellMatCache, ::Any) = ones(1)
        FerriteAssembly.create_cell_state(::CellMatCache, args...) = ones(2)
        grid = generate_grid(Quadrilateral, (10,10))
        addcellset!(grid, "left", x -> x[1]<0.0)
        addcellset!(grid, "right", setdiff(1:getncells(grid), getcellset(grid, "left")))
        ip = Lagrange{RefQuadrilateral,1}()
        dh = DofHandler(grid); add!(dh, :u, ip); close!(dh);
        cv = CellValues(QuadratureRule{RefQuadrilateral}(2), ip)
        userdata = [1.0]
        material = CellMatCache()
        domains = Dict(key => DomainSpec(dh, material, cv; set=getcellset(grid, key), user_data=userdata) for key in ("left", "right"))
        buffers = setup_domainbuffers(domains)
        buffers_ad = setup_domainbuffers(domains; autodiffbuffer=true)
        aold_value = rand()
        aold = ones(ndofs(dh))*aold_value
        _getdomain(dbs::Dict, key::String) = dbs[key]
        _getdomain(sim::Simulation, key) = Simulation(sim.db[key], sim.a, sim.aold)
        for container in (buffers, buffers_ad, Simulation(buffers, nothing, aold), Simulation(buffers_ad, nothing, aold))
            # Basic access functions
            @test FerriteAssembly.get_dofhandler(container) === dh
            @test FerriteAssembly.get_material(container, "left") === material
            @test FerriteAssembly.get_grid(container) === grid
            Δt = rand()
            FerriteAssembly.set_time_increment!(container, Δt)
            # Check properties of stored values
            cont1 = _getdomain(container, "left")
            cellset = FerriteAssembly.getset(cont1)
            @test FerriteAssembly.getset(container, "left") === cellset
            cell_id = first(cellset)
            cb1 = FerriteAssembly.get_itembuffer(cont1)
            sim = isa(container, Simulation) ? cont1 : Simulation(cont1, nothing, aold)
            FerriteAssembly.reinit_buffer!(cb1, sim, cell_id)
            @test FerriteAssembly.get_user_data(cb1) === userdata 
            @test FerriteAssembly.get_user_cache(cb1) == [1.0]
            ae_old = FerriteAssembly.get_aeold(cb1)
            @test all(aold_value == aeold for aeold in ae_old)
            @test cell_id == Ferrite.cellid(cb1)
            @test FerriteAssembly.get_old_state(cb1) === FerriteAssembly.get_old_state(cont1, cell_id)
            @test FerriteAssembly.get_state(cb1) === FerriteAssembly.get_state(cont1, cell_id)
            @test FerriteAssembly.get_old_state(container, "left") === FerriteAssembly.get_old_state(cont1)
            @test FerriteAssembly.get_state(container, "left") === FerriteAssembly.get_state(cont1)
            @test FerriteAssembly.get_state(container)["left"] === FerriteAssembly.get_state(cont1)
            @test FerriteAssembly.get_time_increment(cb1) == Δt

            cont2 = _getdomain(container, "right")
            cb2 = FerriteAssembly.get_itembuffer(cont2)
            @test FerriteAssembly.get_user_cache(cb1) !== FerriteAssembly.get_user_cache(cb2)
            @test FerriteAssembly.get_time_increment(cb2) == Δt
        end
        # Warnings for user input mistakes 
        # Missing cells (hit fewer cells + cells not included)
        @test_logs (:warn,) (:warn,) setup_domainbuffers(Dict("1" => DomainSpec(dh, CellMatCache(), cv; set = getcellset(grid, "left"))))
        # Overlapping sets (hit cells not included)
        @test_logs (:warn,) setup_domainbuffers(Dict(
            "1" => DomainSpec(dh, CellMatCache(), cv; set = Set(1:50)),
            "2" => DomainSpec(dh, CellMatCache(), cv; set = Set(1:50))))
        # Repeated sets (hit more cells)
        @test_logs (:warn,) setup_domainbuffers(Dict(
            "1" => DomainSpec(dh, CellMatCache(), cv),
            "2" => DomainSpec(dh, CellMatCache(), cv)))
    end
    
    @testset "facets" begin
        struct FacetMatCache end 
        FerriteAssembly.allocate_facet_cache(::FacetMatCache, ::Any) = ones(1)
        grid = generate_grid(Quadrilateral, (10,10), Vec((0.0, 0.0)), Vec((1.0, 1.0)))
        ip = Lagrange{RefQuadrilateral,1}()
        dh = DofHandler(grid); add!(dh, :u, ip); close!(dh);
        fv = FacetValues(FacetQuadratureRule{RefQuadrilateral}(2), ip)
        userdata = [1.0]
        domains = Dict(key => DomainSpec(dh, FacetMatCache(), fv; set=getfacetset(grid, key), user_data=userdata) for key in ("left", "right", "top", "bottom"))
        buffers = setup_domainbuffers(domains)
        buffer = buffers["left"]
        aold_value = rand()
        aold = ones(ndofs(dh))*aold_value
        facetbuffer = FerriteAssembly.get_itembuffer(buffer)
        facet_id = first(FerriteAssembly.getset(buffer))
        FerriteAssembly.reinit_buffer!(facetbuffer, Simulation(buffer, zeros(ndofs(dh)), aold), facet_id)
        @test FerriteAssembly.get_user_data(facetbuffer) === userdata 
        @test FerriteAssembly.get_user_cache(facetbuffer) == [1.0]
        @test FerriteAssembly.get_user_cache(facetbuffer) !== FerriteAssembly.get_user_cache(FerriteAssembly.get_itembuffer(buffers["right"]))
        ae_old = FerriteAssembly.get_aeold(facetbuffer)
        @test all(aold_value == aeold for aeold in ae_old)
        @test facet_id[1] == Ferrite.cellid(facetbuffer)
        
        # Not implemented
        @test_throws "AutoDiffBuffer not implemented for FacetBuffer" setup_domainbuffer(domains["left"]; autodiffbuffer=true)

    end

    @testset "mixed cell and facet domains rejected (BUG-011)" begin
        grid = generate_grid(Quadrilateral, (4,4))
        ip = Lagrange{RefQuadrilateral,1}()
        dh = DofHandler(grid); add!(dh, :u, ip); close!(dh)
        cv = CellValues(QuadratureRule{RefQuadrilateral}(2), ip)
        fv = FacetValues(FacetQuadratureRule{RefQuadrilateral}(2), ip)
        cell_domain = DomainSpec(dh, nothing, cv)
        facet_domain = DomainSpec(dh, nothing, fv; set=getfacetset(grid, "left"))

        # Rejection must not depend on which domain type happens to be seen first, so try
        # several dictionary key arrangements (Dict iteration order depends on key hashes).
        key_pairs = (("cells", "facets"), ("facets", "cells"), ("a", "b"), ("z", "y"))
        for (cell_key, facet_key) in key_pairs
            dbs_spec = Dict(cell_key => cell_domain, facet_key => facet_domain)
            @test_throws ArgumentError setup_domainbuffers(dbs_spec)
        end

        # Combining cell and facet integrals instead requires two separate dictionaries.
        cell_dbs = setup_domainbuffers(Dict("cells" => cell_domain))
        facet_dbs = setup_domainbuffers(Dict("facets" => facet_domain))
        v = SimpleIntegrator(Returns(1.0), 0.0)
        work!(v, cell_dbs)
        work!(v, facet_dbs)
        v_cell = SimpleIntegrator(Returns(1.0), 0.0)
        work!(v_cell, cell_dbs)
        v_facet = SimpleIntegrator(Returns(1.0), 0.0)
        work!(v_facet, facet_dbs)
        @test v.val ≈ v_cell.val + v_facet.val
    end

    @testset "empty domains (BUG-016)" begin
        grid = generate_grid(Quadrilateral, (4,4))
        ip = Lagrange{RefQuadrilateral,1}()
        dh = DofHandler(grid); add!(dh, :u, ip); close!(dh)
        cv = CellValues(QuadratureRule{RefQuadrilateral}(2), ip)

        # Explicitly empty set
        @test_throws ArgumentError setup_domainbuffer(DomainSpec(dh, nothing, cv; set=Int[]))
        # Nonempty set whose intersection with the SubDofHandler's cellset is empty
        @test_throws ArgumentError setup_domainbuffer(DomainSpec(dh, nothing, cv; set=[getncells(grid)+1, getncells(grid)+2]))

        # An entirely empty domain dictionary is a no-op, consistent with BUG-008's
        # fix for empty domain dictionaries in work!/load-handler dispatch.
        dbs = setup_domainbuffers(Dict{String, FerriteAssembly.DomainSpec}())
        @test isempty(dbs)
    end
end