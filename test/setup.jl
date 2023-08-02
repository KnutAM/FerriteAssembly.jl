@testset "setup" begin
    @testset "faces" begin
        struct FaceMatCache end 
        FerriteAssembly.allocate_face_cache(::FaceMatCache, ::Any) = ones(1)
        grid = generate_grid(Quadrilateral, (10,10), Vec((0.0, 0.0)), Vec((1.0, 1.0)))
        dh = DofHandler(grid); add!(dh, :u, 1); close!(dh);
        fv = FaceScalarValues(QuadratureRule{1,RefCube}(2), Lagrange{2,RefCube,1}())
        userdata = [1.0]
        domains = Dict(key => DomainSpec(dh, FaceMatCache(), fv; set=getfaceset(grid, key), user_data=userdata) for key in ("left", "right", "top", "bottom"))
        buffers = setup_domainbuffers(domains)
        buffer = buffers["left"]
        aold_value = rand()
        aold = ones(ndofs(dh))*aold_value
        facebuffer = FerriteAssembly.get_itembuffer(buffer)
        face_id = first(FerriteAssembly.getset(buffer))
        FerriteAssembly.reinit_buffer!(facebuffer, buffer, face_id; a=zeros(ndofs(dh)), aold=aold)
        @test FerriteAssembly.get_user_data(facebuffer) === userdata 
        @test FerriteAssembly.get_user_cache(facebuffer) == [1.0]
        @test FerriteAssembly.get_user_cache(facebuffer) !== FerriteAssembly.get_user_cache(FerriteAssembly.get_itembuffer(buffers["right"]))
        ae_old = FerriteAssembly.get_aeold(facebuffer)
        @test all(aold_value == aeold for aeold in ae_old)
        @test face_id[1] == Ferrite.cellid(facebuffer)
        
        # Not implemented
        @test_throws "AutoDiffBuffer not implemented for FaceBuffer" setup_domainbuffer(domains["left"]; autodiffbuffer=true)

    end
end