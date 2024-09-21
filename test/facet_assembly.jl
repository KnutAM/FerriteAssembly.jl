@testset "FacetResidual" begin
    nx, ny = (10, 10)
    grid = generate_grid(Quadrilateral, (nx, ny))
    qr = FacetQuadratureRule{RefQuadrilateral}(2)
    ip = Lagrange{RefQuadrilateral,1}()
    ipg = Lagrange{RefQuadrilateral,1}()
    dh = DofHandler(grid); add!(dh, :u, ip); close!(dh)
    fv = FacetValues(qr, ip, ipg)
    a = zeros(ndofs(dh))
    # Add values only at the boundary
    apply_analytical!(a, dh, :u, x-> x[1] ≈ 1.0 ? 1.0 : 0.0)
    struct FacetMaterial end 
    function FerriteAssembly.facet_residual!(re, ae, ::FacetMaterial, fv, fb)
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
    set = getfacetset(grid, "right")
    buffer = setup_domainbuffer(DomainSpec(dh, FacetMaterial(), fv; set=set))
    r = zeros(ndofs(dh))
    worker = ReAssembler(r)
    work!(worker, buffer; a=a)
    @test sum(r) ≈ 2*sum(a)/(length(set)+1)

    # Full assembler 
    function FerriteAssembly.facet_routine!(Ke, re, ae, ::FacetMaterial, fv::AbstractFacetValues, facetbuffer)
        dofrange = dof_range(facetbuffer, :u)
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
    K = allocate_matrix(dh)
    assembler = start_assemble(K, r2)
    work!(assembler, buffer; a=a2)
    @test K*(a2-a) ≈ (r2-r)
    K3 = allocate_matrix(dh)
    r3 = similar(r)
    kra = KeReAssembler(K3, r3)
    work!(kra, buffer; a=a2)
    @test r3 ≈ r2
    @test K3 ≈ K 
end