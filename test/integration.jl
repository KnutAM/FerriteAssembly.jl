module TestIntegrators
    using FerriteAssembly
    # Re-implement SimpleIntegrator, but with the general Integrator interface. 
    struct MySimpleIntegrand{SI<:SimpleIntegrator}
        integrator::SI
    end
    function FerriteAssembly.integrate_cell!(val::MySimpleIntegrand, cell_state, ae, _, cv, cellbuffer)
        FerriteAssembly.integrate_cell!(val.integrator, cell_state, ae, nothing, cv, cellbuffer)
    end
    function FerriteAssembly.integrate_face!(val::MySimpleIntegrand, ae, _, fv, buffer)
        FerriteAssembly.integrate_face!(val.integrator, ae, nothing, fv, buffer)
    end
end

import .TestIntegrators: MySimpleIntegrand

@testset "Integration" begin
    function solve_problem(dh, ch, buffer; Δt=NaN)
        K = create_sparsity_pattern(dh)
        r = zeros(ndofs(dh))
        a = zeros(ndofs(dh))
        aold = zeros(ndofs(dh))
        assembler = start_assemble(K, r)
        apply!(a, ch)
        set_time_increment!(buffer, Δt)
        work!(assembler, buffer; a=a, aold=aold)
        apply_zero!(K, r, ch)
        norm_r0 = norm(r)
        a .-= K\r 
        apply!(a, ch) # Exact BC
        r_assembler = ReAssembler(r)
        work!(r_assembler, buffer; a=a, aold=aold)
        apply_zero!(r, ch)
        @test norm(r) < norm_r0*1e-12
        return a
    end
    @testset "HeatFlow" begin
        function heatflow_solution(;lx=1.0, ly=1.0, lz=1.0, threading=Val(false))
            grid = generate_grid(Hexahedron, (3,3,3), zero(Vec{3}), Vec(lx, ly, lz))
            dh = DofHandler(grid); add!(dh, :u, 1); close!(dh)
            ch = ConstraintHandler(dh)
            add!(ch, Dirichlet(:u, getfaceset(grid, "left"), Returns(0.0))) 
            add!(ch, Dirichlet(:u, getfaceset(grid, "right"), Returns(1.0)))
            close!(ch)
            update!(ch, 0.0)
            m = EE.StationaryFourier(1.0)
            cv = CellScalarValues(QuadratureRule{3,RefCube}(2), Lagrange{3,RefCube,1}())
            buffer = setup_domainbuffer(DomainSpec(dh, m, cv); threading=threading)
            a = solve_problem(dh, ch, buffer)
            return dh, a, buffer
        end

        function test_heatflow(;threading=Val(false))
            lx = 1.2; ly = 1.3; lz=1.1
            volume = lx*ly*lz
            dh, a, buffer = heatflow_solution(;lx=lx, ly=ly, lz=lz, threading=threading)
            # Calculate the volume 
            volume_integrator = SimpleIntegrator(Returns(1.0), 0.0)
            work!(volume_integrator, buffer)
            @test volume_integrator.val ≈ volume

            # Calculate the average temperature 
            dofsum_integrator = SimpleIntegrator((u, ∇u, state)->u, 0.0)
            work!(dofsum_integrator, buffer; a=a)
            avg_value = dofsum_integrator.val/volume
            @test avg_value ≈ 0.5 # T=0 (left) and T=1 (right)
    
            # Calculate the average temperature gradient
            gradsum_integrator = SimpleIntegrator((u, ∇u, state)->∇u, zero(Vec{3}))
            work!(gradsum_integrator, buffer; a=a)
            avg_grad = gradsum_integrator.val/volume
            @test avg_grad ≈ Vec((1.0/lx, 0.0, 0.0)) # T=0 (left) and T=1 (right)
    
            # Calculate the values above, but together as a heterogeneous tuple
            tuple_integrator = SimpleIntegrator((u, ∇u, state)->(1.0, u, ∇u), (0.0, 0.0, zero(Vec{3})))
            work!(tuple_integrator, buffer; a=a)
            @test volume_integrator.val ≈ tuple_integrator.val[1]
            @test dofsum_integrator.val ≈ tuple_integrator.val[2]
            @test gradsum_integrator.val ≈ tuple_integrator.val[3]
        end
    
        test_heatflow(;threading=Val(false))
        test_heatflow(;threading=Val(true))
        
    end

    @testset "Elasticity" begin
        function elasticity_solution(;lx=1.0, ly=1.0, threading=Val(false))
            grid = generate_grid(Triangle, (3,3), zero(Vec{2}), Vec(lx, ly))
            dh = DofHandler(grid); add!(dh, :u, 2); close!(dh)
            ch = ConstraintHandler(dh)
            add!(ch, Dirichlet(:u, getfaceset(grid, "left"), Returns(zero(Vec{2}))))
            add!(ch, Dirichlet(:u, getfaceset(grid, "right"), Returns(Vec((0.0, 0.1)))))
            close!(ch)
            update!(ch, 0.0)
            m = EE.ElasticPlaneStrain(;E=100.0, ν=0.3)
            cv = CellVectorValues(QuadratureRule{2,RefTetrahedron}(2), Lagrange{2,RefTetrahedron,1}())
            buffer = setup_domainbuffer(DomainSpec(dh, m, cv); threading=threading)
            a = solve_problem(dh, ch, buffer)
            return dh, a, buffer
        end

        function test_elasticity(;threading)
            lx = 1.2; ly = 1.3
            area = lx*ly
            dh, a, buffer = elasticity_solution(;lx=lx, ly=ly, threading=threading)
            # Calculate the area 
            area_integrator = SimpleIntegrator(Returns(1.0), 0.0)
            work!(area_integrator, buffer)
            @test area_integrator.val ≈ area
            # Calculate the average displacements 
            dofsum_integrator = SimpleIntegrator((u, ∇u, state)->u, zero(Vec{2}))
            work!(dofsum_integrator, buffer; a=a)
            avg_value = dofsum_integrator.val/area
            @test abs(avg_value[1]) < 1e-14 # ux=0 (left and right)
            @test avg_value[2] ≈ 0.1/2  # uy=0 (left) and uy=0.1 (right)
    
            # Calculate the average strains gradient
            gradsum_integrator = SimpleIntegrator((u, ∇u, state)->symmetric(∇u), zero(SymmetricTensor{2,2}))
            work!(gradsum_integrator, buffer; a=a)
            avg_grad = gradsum_integrator.val/area
            @test abs(avg_grad[1,1]) < 1e-14
            
            # Calculate the values above, but together as a heterogeneous tuple
            tuple_integrator = SimpleIntegrator((u, ∇u, state)->(1.0, u, symmetric(∇u)), (0.0, zero(Vec{2}), zero(SymmetricTensor{2,2})))
            work!(tuple_integrator, buffer; a=a)
            @test area_integrator.val ≈ tuple_integrator.val[1]
            @test dofsum_integrator.val ≈ tuple_integrator.val[2]
            @test gradsum_integrator.val ≈ tuple_integrator.val[3]
        end
    
        test_elasticity(;threading=Val(false))
        test_elasticity(;threading=Val(true))
    end

    @testset "Poroelasticity" begin
        function poroelasticity_solution(;lx=1.0, ly=1.0, threading=Val(false), ν=0.3)
            ip_u = Lagrange{2,RefTetrahedron,2}()
            ip_p = Lagrange{2,RefTetrahedron,1}()
            ip_geo = Lagrange{2,RefTetrahedron,1}()
            grid = generate_grid(Triangle, (3,3), zero(Vec{2}), Vec(lx, ly))
            dh = DofHandler(grid); add!(dh, :u, 2, ip_u); add!(dh, :p, 1, ip_p); close!(dh)
            ch = ConstraintHandler(dh)
            add!(ch, Dirichlet(:u, getfaceset(grid, "left"), Returns(zero(Vec{2}))))
            add!(ch, Dirichlet(:u, getfaceset(grid, "right"), Returns(Vec((0.1, 0.0)))))
            add!(ch, Dirichlet(:p, getfaceset(grid, "left"), Returns(0.0)))
            add!(ch, Dirichlet(:p, getfaceset(grid, "right"), Returns(1.0)))
            close!(ch)
            update!(ch, 0.0)
            m = EE.PoroElasticPlaneStrain(;E=100.0, ν=ν, k=1.0, 
                α=0.0, # Ensure no coupling between u and p
                β=0.0  # Remove time-dependence on p. Time-dependence on u only via coupling. 
                )
            qr = QuadratureRule{2,RefTetrahedron}(2)
            cv = (u=CellVectorValues(qr, ip_u, ip_geo), p=CellScalarValues(qr, ip_p, ip_geo))
            buffer = setup_domainbuffer(DomainSpec(dh, m, cv); threading=threading)
            a = solve_problem(dh, ch, buffer; Δt=1.0)
            return dh, a, buffer
        end

        function test_poroelasticity(;threading)
            lx = 1.8; ly = 1.7; ν=0.3
            area = lx*ly
            dh, a, buffer = poroelasticity_solution(;lx=lx, ly=ly, threading=threading, ν=ν)
            # Calculate the area 
            area_integrator = SimpleIntegrator(Returns(1.0), 0.0)
            work!(area_integrator, buffer)
            @test area_integrator.val ≈ area
            # Calculate the average displacements for tension with fixed left end
            # I.e. not uniaxial, but symmetric across horizontal center line
            u_integrator = SimpleIntegrator((u, ∇u, state)->u[:u], zero(Vec{2}))
            work!(u_integrator, buffer; a=a)
            u_avg = u_integrator.val/area
            @test u_avg[1] ≈ 0.1/2      # ux=0 (left) and ux=0.1 (right)
            @test abs(u_avg[2]) < 1e-14 # uy=0 (left) and uy "free" (right) (symmetric)

            ∇u_integrator = SimpleIntegrator((u, ∇u, state)->symmetric(∇u[:u]), zero(SymmetricTensor{2,2}))
            work!(∇u_integrator, buffer; a=a)
            ∇u_avg = ∇u_integrator.val/area
            ϵxx_avg = 0.1/lx            # ux=0 (left) and ux=0.1 (right)
            @test ∇u_avg[1,1] ≈ ϵxx_avg
            ϵyy_avg = -ϵxx_avg*ν # uy=0 (left) and uy "free" (right) (symmetric)
            @test isapprox(∇u_avg[2,2], ϵyy_avg; rtol=0.1) # Only approximately
            @test abs(∇u_avg[1,2]) < 1e-2 # Should cancel due to symmetry, but only approximately
            
            # Calculate the average pressure
            p_integrator = SimpleIntegrator((u, ∇u, state)->u[:p], 0.0)
            work!(p_integrator, buffer; a=a)
            p_avg = p_integrator.val/area
            @test p_avg ≈ 0.5 # p=0 (left) and p=1 (right)
    
            # Calculate the average pressure gradient
            ∇p_integrator = SimpleIntegrator((u, ∇u, state)->∇u[:p], zero(Vec{2}))
            work!(∇p_integrator, buffer; a=a)
            ∇p_avg = ∇p_integrator.val/area
            @test ∇p_avg ≈ Vec((1.0/lx, 0.0)) # T=0 (left) and T=1 (right)
            
            # Calculate values from above, but together as a heterogeneous tuple
            tuple_integrator = SimpleIntegrator((u, ∇u, state)->(1.0, u[:p], symmetric(∇u[:u])), (0.0, 0.0, zero(SymmetricTensor{2,2})))
            work!(tuple_integrator, buffer; a=a)
            @test area_integrator.val ≈ tuple_integrator.val[1]
            @test p_integrator.val ≈ tuple_integrator.val[2]
            @test ∇u_integrator.val ≈ tuple_integrator.val[3]

            # Do the same test, but using Integrator with MySimpleIntegrand
            simple_integrator = SimpleIntegrator((u, ∇u, state)->(1.0, u[:p], symmetric(∇u[:u])), (0.0, 0.0, zero(SymmetricTensor{2,2})))
            integrator = Integrator(MySimpleIntegrand(simple_integrator))
            work!(integrator, buffer; a=a)
            @test all(simple_integrator.val .≈ tuple_integrator.val)
        end

        test_poroelasticity(;threading=Val(false))
        test_poroelasticity(;threading=Val(true)) 

    end

    @testset "MultipleDomains" begin
        function multidomain_heatflow_solution(;lx=1.0, ly=1.0, lz=1.0, threading=Val(false))
            grid = generate_grid(Hexahedron, (3,3,3), zero(Vec{3}), Vec(lx, ly, lz))
            addcellset!(grid, "set1", x->x[1]<lx/2)
            addcellset!(grid, "set2", setdiff!(Set(1:getncells(grid)), getcellset(grid, "set1")))

            dh = DofHandler(grid); add!(dh, :u, 1); close!(dh)
            ch = ConstraintHandler(dh)
            add!(ch, Dirichlet(:u, getfaceset(grid, "left"), Returns(0.0))) 
            add!(ch, Dirichlet(:u, getfaceset(grid, "right"), Returns(1.0)))
            close!(ch)
            update!(ch, 0.0)
            m = EE.StationaryFourier(1.0)
            cv = CellScalarValues(QuadratureRule{3,RefCube}(2), Lagrange{3,RefCube,1}())
            ad1 = DomainSpec(dh, m, cv; set=getcellset(grid, "set1"))
            ad2 = DomainSpec(dh, m, cv; set=getcellset(grid, "set2"))
            buffers = setup_domainbuffers(Dict("d1"=>ad1, "d2"=>ad2); threading=threading)
            a = solve_problem(dh, ch, buffers)
            return dh, a, buffers
        end

        function test_multidomain_heatflow(;threading=Val(false))
            lx = 1.2; ly = 1.3; lz=1.1
            volume = lx*ly*lz
            dh, a, buffer = multidomain_heatflow_solution(;lx=lx, ly=ly, lz=lz, threading=threading)
            volume1 = volume*length(getcellset(dh.grid, "set1"))/getncells(dh.grid)
            volume2 = volume*length(getcellset(dh.grid, "set2"))/getncells(dh.grid)
            
            # Calculate the full volume with SimpleIntegrator
            volume_integrator = SimpleIntegrator(Returns(1.0), 0.0)
            work!(volume_integrator, buffer)
            @test volume_integrator.val ≈ volume

            # Calculate the volumes for each domain 
            v1_integrator = SimpleIntegrator(Returns(1.0), 0.0; domains="d1")
            work!(v1_integrator, buffer)
            v1 = v1_integrator.val
            @test v1 ≈ volume1
            v2_integrator = SimpleIntegrator(Returns(1.0), 0.0; domains="d2")
            work!(v2_integrator, buffer)
            v2 = v2_integrator.val
            @test v2 ≈ volume2
            @test (v1 + v2) ≈ volume 
            
            # Calculate the volumes for both domains by specifying them
            v_integrator = SimpleIntegrator(Returns(1.0), 0.0; domains=("d1", "d2"))
            work!(v_integrator, buffer)
            @test v_integrator.val ≈ volume 

            # Check v1 using Integrator (with MySimpleIntegrand that wraps a SimpleIntegrator)
            v1_fullintegrator = Integrator(MySimpleIntegrand(SimpleIntegrator(Returns(1.0), 0.0)); domains="d1")
            work!(v1_fullintegrator, buffer)
            @test v1_fullintegrator.val.integrator.val ≈ v1
        end
    
        test_multidomain_heatflow(;threading=Val(false))
        test_multidomain_heatflow(;threading=Val(true))
    end

end

@testset "FaceIntegration" begin
    @testset "Geometric tests" begin
        w, h = rand(2)
        grid = generate_grid(Quadrilateral, (10,10), Vec((0.0, 0.0)), Vec((w, h)))
        dh = DofHandler(grid); add!(dh, :u, 1); close!(dh);
        fv = FaceScalarValues(QuadratureRule{1,RefCube}(2), Lagrange{2,RefCube,1}())
        domains = Dict(key => DomainSpec(dh, nothing, fv; set=getfaceset(grid, key)) for key in ("left", "right", "top", "bottom"))
        buffers = setup_domainbuffers(domains)
        for (Aref, domain_keys) in ((2*(w+h), nothing), (2*w, ("top", "bottom"), (h, "left")))
            ig = SimpleIntegrator((u,∇u,n)->1.0, 0.0; domains=domain_keys)
            work!(ig, buffers)
            @test ig.val ≈ Aref
            ig.val = 0.0
            full_ig = Integrator(MySimpleIntegrand(ig); domains=domain_keys)
            work!(full_ig, buffers)
            @test ig.val ≈ Aref
        end
    @testset "Function values" begin
        lx, ly, lz = rand(3)
        p = Vec((lx, ly, lz))/2
        grid = generate_grid(Tetrahedron, (10,10,10), -p, p)
        dh = DofHandler(grid); add!(dh, :u, 3); add!(dh, :v, 1); close!(dh);
        qr = QuadratureRule{2,RefTetrahedron}(2)
        ip = Lagrange{3,RefTetrahedron,1}()
        fv_u = FaceVectorValues(qr, ip)
        fv_v = FaceScalarValues(qr, ip)

        a = zeros(ndofs(dh))
        apply_analytical!(a, dh, :u, identity) # set dofs equal to coordinates
        apply_analytical!(a, dh, :v, first)  # 
        domains = Dict(key => DomainSpec(dh, nothing, (u=fv_u, v=fv_v); set=getfaceset(grid, key)) for key in ("left", "right", "top", "bottom"))
        buffers = setup_domainbuffers(domains)

        ig = SimpleIntegrator((u,∇u,n)->(1.0, u.u ⋅ n, u.v), (0.0, 0.0, 0.0))
        work!(ig, buffers; a=a)
        @test ig.val[1] ≈ ly*2*(lx+lz)
        @test ig.val[2] ≈ 2*lx*ly*lz
        @test isapprox(ig.val[3], 0.0; atol=ig.val[1]*1e-8)
    end
    end
end