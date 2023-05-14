module TestIntegrators
    using FerriteAssembly
    # Re-implement SimpleIntegrator, but with the general Integrator interface. 
    struct MySimpleIntegrand{SI<:SimpleIntegrator}
        integrator::SI
    end
    function FerriteAssembly.integrate_cell!(val::MySimpleIntegrand, cell_state, ae, _, cv, cellbuffer)
        FerriteAssembly.simple_integrate_cell!(val.integrator, cell_state, ae, cv, cellbuffer)
    end
end

import .TestIntegrators: MySimpleIntegrand

@testset "Integration" begin
    function solve_problem(dh, ch, buffer, states; Δt=NaN)
        K = create_sparsity_pattern(dh)
        r = zeros(ndofs(dh))
        a = zeros(ndofs(dh))
        aold = zeros(ndofs(dh))
        assembler = start_assemble(K, r)
        apply!(a, ch)
        doassemble!(assembler, states, buffer; a=a, aold=aold, Δt=Δt)
        apply_zero!(K, r, ch)
        norm_r0 = norm(r)
        a .-= K\r 
        apply!(a, ch) # Exact BC
        r_assembler = ReAssembler(r)
        doassemble!(r_assembler, states, buffer; a=a, aold=aold, Δt=Δt)
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
            buffer, states = setup_assembly(dh, m, cv; threading=threading)
            a = solve_problem(dh, ch, buffer, states)
            return dh, a, states, buffer
        end

        function test_heatflow(;threading=Val(false))
            lx = 1.2; ly = 1.3; lz=1.1
            volume = lx*ly*lz
            dh, a, states, buffer = heatflow_solution(;lx=lx, ly=ly, lz=lz, threading=threading)
            # Calculate the volume 
            volume_integrator = SimpleIntegrator(Returns(1.0), 0.0)
            doassemble!(volume_integrator, states, buffer)
            @test volume_integrator.val ≈ volume

            # Calculate the average temperature 
            dofsum_integrator = SimpleIntegrator((u, ∇u, state)->u, 0.0)
            doassemble!(dofsum_integrator, states, buffer; a=a)
            avg_value = dofsum_integrator.val/volume
            @test avg_value ≈ 0.5 # T=0 (left) and T=1 (right)
    
            # Calculate the average temperature gradient
            gradsum_integrator = SimpleIntegrator((u, ∇u, state)->∇u, zero(Vec{3}))
            doassemble!(gradsum_integrator, states, buffer; a=a)
            avg_grad = gradsum_integrator.val/volume
            @test avg_grad ≈ Vec((1.0/lx, 0.0, 0.0)) # T=0 (left) and T=1 (right)
    
            # Calculate the values above, but together as a heterogeneous tuple
            tuple_integrator = SimpleIntegrator((u, ∇u, state)->(1.0, u, ∇u), (0.0, 0.0, zero(Vec{3})))
            doassemble!(tuple_integrator, states, buffer; a=a)
            @test volume_integrator.val == tuple_integrator.val[1]
            @test dofsum_integrator.val == tuple_integrator.val[2]
            @test gradsum_integrator.val == tuple_integrator.val[3]
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
            buffer, states = setup_assembly(dh, m, cv; threading=threading)
            a = solve_problem(dh, ch, buffer, states)
            return dh, a, states, buffer
        end

        function test_elasticity(;threading)
            lx = 1.2; ly = 1.3
            area = lx*ly
            dh, a, states, buffer = elasticity_solution(;lx=lx, ly=ly, threading=threading)
            # Calculate the area 
            area_integrator = SimpleIntegrator(Returns(1.0), 0.0)
            doassemble!(area_integrator, states, buffer)
            @test area_integrator.val ≈ area
            # Calculate the average displacements 
            dofsum_integrator = SimpleIntegrator((u, ∇u, state)->u, zero(Vec{2}))
            doassemble!(dofsum_integrator, states, buffer; a=a)
            avg_value = dofsum_integrator.val/area
            @test abs(avg_value[1]) < 1e-14 # ux=0 (left and right)
            @test avg_value[2] ≈ 0.1/2  # uy=0 (left) and uy=0.1 (right)
    
            # Calculate the average strains gradient
            gradsum_integrator = SimpleIntegrator((u, ∇u, state)->symmetric(∇u), zero(SymmetricTensor{2,2}))
            doassemble!(gradsum_integrator, states, buffer; a=a)
            avg_grad = gradsum_integrator.val/area
            @test abs(avg_grad[1,1]) < 1e-14
            
            # Calculate the values above, but together as a heterogeneous tuple
            tuple_integrator = SimpleIntegrator((u, ∇u, state)->(1.0, u, symmetric(∇u)), (0.0, zero(Vec{2}), zero(SymmetricTensor{2,2})))
            doassemble!(tuple_integrator, states, buffer; a=a)
            @test area_integrator.val == tuple_integrator.val[1]
            @test dofsum_integrator.val == tuple_integrator.val[2]
            @test gradsum_integrator.val == tuple_integrator.val[3]
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
            buffer, states = setup_assembly(dh, m, cv; threading=threading)
            a = solve_problem(dh, ch, buffer, states; Δt=1.0)
            return dh, a, states, buffer
        end

        function test_poroelasticity(;threading)
            lx = 1.8; ly = 1.7; ν=0.3
            area = lx*ly
            dh, a, states, buffer = poroelasticity_solution(;lx=lx, ly=ly, threading=threading, ν=ν)
            # Calculate the area 
            area_integrator = SimpleIntegrator(Returns(1.0), 0.0)
            doassemble!(area_integrator, states, buffer)
            @test area_integrator.val ≈ area
            # Calculate the average displacements for tension with fixed left end
            # I.e. not uniaxial, but symmetric across horizontal center line
            u_integrator = SimpleIntegrator((u, ∇u, state)->u[:u], zero(Vec{2}))
            doassemble!(u_integrator, states, buffer; a=a)
            u_avg = u_integrator.val/area
            @test u_avg[1] ≈ 0.1/2      # ux=0 (left) and ux=0.1 (right)
            @test abs(u_avg[2]) < 1e-14 # uy=0 (left) and uy "free" (right) (symmetric)

            ∇u_integrator = SimpleIntegrator((u, ∇u, state)->symmetric(∇u[:u]), zero(SymmetricTensor{2,2}))
            doassemble!(∇u_integrator, states, buffer; a=a)
            ∇u_avg = ∇u_integrator.val/area
            ϵxx_avg = 0.1/lx            # ux=0 (left) and ux=0.1 (right)
            @test ∇u_avg[1,1] ≈ ϵxx_avg
            ϵyy_avg = -ϵxx_avg*ν # uy=0 (left) and uy "free" (right) (symmetric)
            @test isapprox(∇u_avg[2,2], ϵyy_avg; rtol=0.1) # Only approximately
            @test abs(∇u_avg[1,2]) < 1e-2 # Should cancel due to symmetry, but only approximately
            
            # Calculate the average pressure
            p_integrator = SimpleIntegrator((u, ∇u, state)->u[:p], 0.0)
            doassemble!(p_integrator, states, buffer; a=a)
            p_avg = p_integrator.val/area
            @test p_avg ≈ 0.5 # p=0 (left) and p=1 (right)
    
            # Calculate the average pressure gradient
            ∇p_integrator = SimpleIntegrator((u, ∇u, state)->∇u[:p], zero(Vec{2}))
            doassemble!(∇p_integrator, states, buffer; a=a)
            ∇p_avg = ∇p_integrator.val/area
            @test ∇p_avg ≈ Vec((1.0/lx, 0.0)) # T=0 (left) and T=1 (right)
            
            # Calculate values from above, but together as a heterogeneous tuple
            tuple_integrator = SimpleIntegrator((u, ∇u, state)->(1.0, u[:p], symmetric(∇u[:u])), (0.0, 0.0, zero(SymmetricTensor{2,2})))
            doassemble!(tuple_integrator, states, buffer; a=a)
            @test area_integrator.val == tuple_integrator.val[1]
            @test p_integrator.val == tuple_integrator.val[2]
            @test ∇u_integrator.val == tuple_integrator.val[3]

            # Do the same test, but using Integrator with MySimpleIntegrand
            simple_integrator = SimpleIntegrator((u, ∇u, state)->(1.0, u[:p], symmetric(∇u[:u])), (0.0, 0.0, zero(SymmetricTensor{2,2})))
            integrator = Integrator(MySimpleIntegrand(simple_integrator))
            doassemble!(integrator, states, buffer; a=a)
            if isa(threading, Val{true})
                @test all(simple_integrator.val .≈ tuple_integrator.val) # Integrator doesn't thread, so only inexact equality
            else
                @test simple_integrator.val == tuple_integrator.val # Should be exact 
            end
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
            ad1 = AssemblyDomain("d1", dh, m, cv; cellset=getcellset(grid, "set1"))
            ad2 = AssemblyDomain("d2", dh, m, cv; cellset=getcellset(grid, "set2"))
            buffer, states = setup_assembly([ad1, ad2]; threading=threading)
            a = solve_problem(dh, ch, buffer, states)
            return dh, a, states, buffer
        end

        function test_multidomain_heatflow(;threading=Val(false))
            lx = 1.2; ly = 1.3; lz=1.1
            volume = lx*ly*lz
            dh, a, states, buffer = multidomain_heatflow_solution(;lx=lx, ly=ly, lz=lz, threading=threading)
            volume1 = volume*length(getcellset(dh.grid, "set1"))/getncells(dh.grid)
            volume2 = volume*length(getcellset(dh.grid, "set2"))/getncells(dh.grid)
            
            # Calculate the full volume with SimpleIntegrator
            volume_integrator = SimpleIntegrator(Returns(1.0), 0.0)
            doassemble!(volume_integrator, states, buffer)
            @test volume_integrator.val ≈ volume

            # Calculate the volumes for each domain 
            v1_integrator = SimpleIntegrator(Returns(1.0), 0.0; domains="d1")
            doassemble!(v1_integrator, states, buffer)
            v1 = v1_integrator.val
            @test v1 ≈ volume1
            v2_integrator = SimpleIntegrator(Returns(1.0), 0.0; domains="d2")
            doassemble!(v2_integrator, states, buffer)
            v2 = v2_integrator.val
            @test v2 ≈ volume2
            @test (v1 + v2) ≈ volume 
            
            # Calculate the volumes for both domains by specifying them
            v_integrator = SimpleIntegrator(Returns(1.0), 0.0; domains=("d1", "d2"))
            doassemble!(v_integrator, states, buffer)
            @test v_integrator.val ≈ volume 

            # Check v1 using Integrator (with MySimpleIntegrand that wraps a SimpleIntegrator)
            v1_fullintegrator = Integrator(MySimpleIntegrand(SimpleIntegrator(Returns(1.0), 0.0)); domains="d1")
            doassemble!(v1_fullintegrator, states, buffer)
            @test v1_fullintegrator.val.integrator.val ≈ v1
        end
    
        test_multidomain_heatflow(;threading=Val(false))
        test_multidomain_heatflow(;threading=Val(true))
    end

end