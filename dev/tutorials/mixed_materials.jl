using Ferrite, FerriteAssembly, MaterialModelsBase, WriteVTK
import FerriteAssembly.ExampleElements: J2Plasticity, ElasticPlaneStrain

function create_grid_with_inclusion()
    p1 = Vec((-1.0, -1.0))
    p2 = Vec(( 1.0,  1.0))
    grid = generate_grid(Quadrilateral, (20,20), p1, p2)
    addcellset!(grid, "inclusion", x -> norm(x) < 0.8)
    addcellset!(grid, "matrix", setdiff(1:getncells(grid), getcellset(grid, "inclusion")))
    return grid
end
grid = create_grid_with_inclusion();

ip = Lagrange{RefQuadrilateral,1}()^2;

dh = DofHandler(grid)
add!(dh, :u, ip)
close!(dh);

ch = ConstraintHandler(dh)
add!(ch, Dirichlet(:u, getfacetset(grid,"left"), Returns(zero(Vec{2}))))
f_dbc(x,t) = Vec((0.05*t, 0.0))
add!(ch, Dirichlet(:u, getfacetset(grid, "right"), f_dbc))
close!(ch);

qr = QuadratureRule{RefQuadrilateral}(2)
cv = CellValues(qr, ip);

elastic_material = ElasticPlaneStrain(;E=210e3, ν=0.3)
plastic_material = ReducedStressState(
    PlaneStrain(),
    J2Plasticity(;E=210e3, ν=0.3, σ0=100.0, H=10e3));

domains = Dict(
    "elastic"=>DomainSpec(dh, elastic_material, cv; set=getcellset(grid, "inclusion")),
    "plastic"=>DomainSpec(dh, plastic_material, cv; set=getcellset(grid, "matrix")) );

buffer = setup_domainbuffers(domains);

function solve_nonlinear_timehistory(buffer, dh, ch; time_history)
    maxiter = 10
    tolerance = 1e-6
    K = allocate_matrix(dh)
    r = zeros(ndofs(dh))
    a = zeros(ndofs(dh))
    # Prepare postprocessing
    pvd = paraview_collection("multiple_materials")
    VTKGridFile("multiple_materials_0", dh) do vtk
        write_solution(vtk, dh, a)
        Ferrite.write_cellset(vtk, dh.grid, "inclusion")
        pvd[0.0] = vtk
    end
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
            a .-= K\r
            apply!(a, ch)
        end
        # Postprocess
        VTKGridFile("multiple_materials_$n", dh) do vtk
            write_solution(vtk, dh, a)
            Ferrite.write_cellset(vtk, dh.grid, "inclusion")
            pvd[t] = vtk
        end
        # If converged, update the old state variables to the current.
        update_states!(buffer)
    end
    close(pvd)
    return nothing
end;
solve_nonlinear_timehistory(buffer, dh, ch; time_history=collect(range(0,1,20)));

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
