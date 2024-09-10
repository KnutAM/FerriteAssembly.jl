using Ferrite, FerriteAssembly, MaterialModelsBase, MechanicalMaterialModels, WriteVTK

function create_grid_with_inclusion()
    p1 = Vec((-1.0, -1.0))
    p2 = Vec(( 1.0,  1.0))
    grid = generate_grid(Quadrilateral, (20, 20), p1, p2)
    addcellset!(grid, "inclusion", x -> norm(x) < 0.5)
    addcellset!(grid, "matrix", setdiff(1:getncells(grid), getcellset(grid, "inclusion")))
    return grid
end
grid = create_grid_with_inclusion();

ip = Lagrange{RefQuadrilateral,1}()^2;

dh = DofHandler(grid)
add!(dh, :u, ip)
close!(dh);

ch = ConstraintHandler(dh)
add!(ch, Dirichlet(:u, getfacetset(grid,"left"), Returns(0.0), 1))
add!(ch, Dirichlet(:u, getfacetset(grid,"bottom"), Returns(0.0), 2))
f_dbc(x,t) = 0.02 * t # 1 % strain at t = 1
add!(ch, Dirichlet(:u, getfacetset(grid, "right"), f_dbc, 1))
close!(ch);

qr = QuadratureRule{RefQuadrilateral}(2)
cv = CellValues(qr, ip);

elastic_material = ReducedStressState(PlaneStrain(), LinearElastic(;E=210e3, ν=0.3))

plastic_material = ReducedStressState(PlaneStrain(), Plastic(;
    elastic   = LinearElastic(E = 210e3, ν = 0.3),
    yield     = 100.0,
    isotropic = Voce(;Hiso = 10e3, κ∞ = Inf),
    kinematic = ArmstrongFrederick(;Hkin = 0.0, β∞ = 1.0) # No kinematic hardening
    ));

domains = Dict(
    "elastic"=>DomainSpec(dh, elastic_material, cv; set=getcellset(grid, "inclusion")),
    "plastic"=>DomainSpec(dh, plastic_material, cv; set=getcellset(grid, "matrix")) );

buffer = setup_domainbuffers(domains);

function calculate_stress(m::ReducedStressState, u, ∇u, qp_state)
    ϵ = MaterialModelsBase.expand_tensordim(m.stress_state, symmetric(∇u))
    σ = calculate_stress(m.material, ϵ, qp_state)
    return MaterialModelsBase.reduce_tensordim(m.stress_state, σ)
end
calculate_stress(m::LinearElastic, ϵ, qp_state) = m.C ⊡ ϵ
calculate_stress(m::Plastic, ϵ, qp_state) = calculate_stress(m.elastic, ϵ - qp_state.ϵp, qp_state)

qe = QuadratureEvaluator{SymmetricTensor{2,2,Float64,3}}(buffer, calculate_stress)

proj = L2Projector(grid)
add!(proj, 1:getncells(grid), ip; qr_rhs = qr)
close!(proj)

function solve_nonlinear_timehistory(buffer, dh, ch, l2_proj, qp_evaluator; time_history)
    maxiter = 10
    tolerance = 1e-6
    K = allocate_matrix(dh)
    r = zeros(ndofs(dh))
    a = zeros(ndofs(dh))
    # Prepare postprocessing
    pvd = paraview_collection("multiple_materials")
    for (n, t) in enumerate(time_history)
        # Update and apply the Dirichlet boundary conditions
        update!(ch, t)
        apply!(a, ch)
        for i in 1:maxiter
            # Assemble the system
            assembler = start_assemble(K, r)
            work!(assembler, buffer; a=a)
            # Apply boundary conditions
            apply_zero!(K, r, ch)
            # Check convergence
            norm(r) < tolerance && break
            i == maxiter && error("Did not converge")
            # Solve the linear system and update the dof vector
            a .-= K \ r
            apply!(a, ch)
        end

        # If converged, update the old state variables to the current.
        update_states!(buffer)

        # Postprocess
        work!(qp_evaluator, buffer; a=a)
        stresses = project(l2_proj, qp_evaluator.data)
        VTKGridFile("multiple_materials_$n", dh) do vtk
            write_solution(vtk, dh, a)
            write_projection(vtk, l2_proj, stresses, "stress")
            Ferrite.write_cellset(vtk, dh.grid, "inclusion")
            pvd[t] = vtk
        end
    end
    close(pvd)
    return nothing
end;
solve_nonlinear_timehistory(buffer, dh, ch, proj, qe; time_history=collect(range(0, 1, 40)));

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
