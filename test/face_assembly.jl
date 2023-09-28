@testset "FaceResidual" begin
    nx, ny = (10, 10)
    grid = generate_grid(Quadrilateral, (nx, ny))
    qr = QuadratureRule{1,RefCube}(2)
    ip = Lagrange{2,RefCube,1}()
    ipg = Lagrange{2,RefCube,1}()
    dh = DofHandler(grid); add!(dh, :u, 1, ip); close!(dh)
    fv = FaceScalarValues(qr, ip, ipg)
    a = zeros(ndofs(dh))
    # Add values only at the boundary
    apply_analytical!(a, dh, :u, x-> x[1] ≈ 1.0 ? 1.0 : 0.0)
    struct FaceMaterial end 
    function FerriteAssembly.face_residual!(re, ae, ::FaceMaterial, fv, fb)
        dofrange = dof_range(fb, :u)
        for q_point in 1:getnquadpoints(fv)
            dΓ = getdetJdV(fv, q_point)
            u = function_value(fv, q_point, ae)
            for (i, I) in pairs(dofrange)
                δN = shape_value(fv, q_point, i)
                re[I] += δN*u*dΓ
            end
        end
    end
    set = getfaceset(grid, "right")
    buffer = setup_domainbuffer(DomainSpec(dh, FaceMaterial(), fv; set=set))
    r = zeros(ndofs(dh))
    worker = ReAssembler(r)
    work!(worker, buffer; a=a)
    @test sum(r) ≈ 2*sum(a)/(length(set)+1)

    # Full assembler 
    function FerriteAssembly.face_routine!(Ke, re, ae, ::FaceMaterial, fv::FaceValues, facebuffer)
        dofrange = dof_range(facebuffer, :u)
        for q_point in 1:getnquadpoints(fv)
            dΓ = getdetJdV(fv, q_point)
            u = function_value(fv, q_point, ae)
            for (i, I) in pairs(dofrange)
                δN = shape_value(fv, q_point, i)
                re[I] += δN*u*dΓ
                for (j, J) in pairs(dofrange)
                    N = shape_value(fv, q_point, j)
                    Ke[I,J] += δN*N*dΓ
                end
            end
        end
    end

    a2 = similar(a)
    r2 = similar(r)
    map!(x -> x≈0 ? zero(x) : rand(), a2, a)
    K = create_sparsity_pattern(dh)
    assembler = start_assemble(K, r2)
    work!(assembler, buffer; a=a2)
    @test K*(a2-a) ≈ (r2-r)
    K3 = create_sparsity_pattern(dh)
    r3 = similar(r)
    kra = KeReAssembler(K3, r3)
    work!(kra, buffer; a=a2)
    @test r3 ≈ r2
    @test K3 ≈ K 
end