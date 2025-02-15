@testset "setup" begin
    
    @testset "cells" begin
        struct CellMatCache end 
        FerriteAssembly.allocate_cell_cache(::CellMatCache, ::Any) = ones(1)
        FerriteAssembly.create_cell_state(::CellMatCache, args...) = ones(2)
        grid = generate_grid(Quadrilateral, (10,10))
        addcellset!(grid, "left", x -> x[1]<0.0)
        addcellset!(grid, "right", setdiff(1:getncells(grid), getcellset(grid, "left")))
        ip = Lagrange{RefQuadrilateral,1}()
        dh = DofHandler(grid); add!(dh, :u, ip); close!(dh);
        cv = CellValues(QuadratureRule{RefQuadrilateral}(2), ip)
        userdata = [1.0]
        domains = Dict(key => DomainSpec(dh, CellMatCache(), cv; set=getcellset(grid, key), user_data=userdata) for key in ("left", "right"))
        buffers = setup_domainbuffers(domains)
        buffers_ad = setup_domainbuffers(domains; autodiffbuffer=true)
        aold_value = rand()
        aold = ones(ndofs(dh))*aold_value
        for domainbuffers in (buffers, buffers_ad)
            db1 = domainbuffers["left"]
            cellset = FerriteAssembly.getset(db1)
            @test FerriteAssembly.getset(domainbuffers, "left") === cellset
            cell_id = first(cellset)            
            cb1 = FerriteAssembly.get_itembuffer(db1)
            FerriteAssembly.reinit_buffer!(cb1, db1, cell_id; a=zeros(ndofs(dh)), aold=aold)
            @test FerriteAssembly.get_user_data(cb1) === userdata 
            @test FerriteAssembly.get_user_cache(cb1) == [1.0]
            ae_old = FerriteAssembly.get_aeold(cb1)
            @test all(aold_value == aeold for aeold in ae_old)
            @test cell_id == Ferrite.cellid(cb1)
            @test FerriteAssembly.get_old_state(cb1) === FerriteAssembly.get_old_state(db1, cell_id)
            @test FerriteAssembly.get_state(cb1) === FerriteAssembly.get_state(db1, cell_id)
            @test FerriteAssembly.get_old_state(domainbuffers, "left") === FerriteAssembly.get_old_state(db1)
            @test FerriteAssembly.get_state(domainbuffers, "left") === FerriteAssembly.get_state(db1)

            db2 = domainbuffers["right"]
            cb2 = FerriteAssembly.get_itembuffer(db2)
            @test FerriteAssembly.get_user_cache(cb1) !== FerriteAssembly.get_user_cache(cb2)
        end
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
        FerriteAssembly.reinit_buffer!(facetbuffer, buffer, facet_id; a=zeros(ndofs(dh)), aold=aold)
        @test FerriteAssembly.get_user_data(facetbuffer) === userdata 
        @test FerriteAssembly.get_user_cache(facetbuffer) == [1.0]
        @test FerriteAssembly.get_user_cache(facetbuffer) !== FerriteAssembly.get_user_cache(FerriteAssembly.get_itembuffer(buffers["right"]))
        ae_old = FerriteAssembly.get_aeold(facetbuffer)
        @test all(aold_value == aeold for aeold in ae_old)
        @test facet_id[1] == Ferrite.cellid(facetbuffer)
        
        # Not implemented
        @test_throws "AutoDiffBuffer not implemented for FacetBuffer" setup_domainbuffer(domains["left"]; autodiffbuffer=true)

    end
end