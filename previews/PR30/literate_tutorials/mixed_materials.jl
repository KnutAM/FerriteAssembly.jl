# # Multiple materials
# In this tutorial, we will see how we can assemble a domain and solve a problem 
# where we have multiple material behaviors. In this simple case, we will consider 
# an elastic inclusion, embedded in a plastically deforming matrix. 
# For this example, we'll use material models defined in the 
# ![results](mixed_materials.png)
# **Figure 1:** Results showing the norm of the stress
# 
# [`MechanicalMaterialModels.jl`](https://github.com/KnutAM/MechanicalMaterialModels.jl)
# package, which defines models according to the `MaterialModelsBase` interface.
# 
# The full script without intermediate comments is available at the 
# [bottom of this page](@ref mixed_materials_plain_program).
#
using Ferrite, FerriteAssembly, MaterialModelsBase, MechanicalMaterialModels, WriteVTK

# ## Setup Ferrite quantities
# We start by the grid with sets for the inclusion with 
# radius 0.5 and the surrounding matrix. 
function create_grid_with_inclusion()
    p1 = Vec((-1.0, -1.0))
    p2 = Vec(( 1.0,  1.0))
    grid = generate_grid(Quadrilateral, (100, 100), p1, p2)
    grid = generate_grid(Quadrilateral, (10, 10), p1, p2)   #src
    addcellset!(grid, "inclusion", x -> norm(x) < 0.5)
    addcellset!(grid, "matrix", setdiff(1:getncells(grid), getcellset(grid, "inclusion")))
    return grid 
end
grid = create_grid_with_inclusion();

# Define interpolation
ip = Lagrange{RefQuadrilateral,1}()^2;

# Followed by the dof handler 
dh = DofHandler(grid)
add!(dh, :u, ip)
close!(dh);

# And then Dirichlet conditions
ch = ConstraintHandler(dh)
add!(ch, Dirichlet(:u, getfacetset(grid,"left"), Returns(0.0), 1))
add!(ch, Dirichlet(:u, getfacetset(grid,"bottom"), Returns(0.0), 2))
f_dbc(x,t) = 0.01 * t # 1 % strain at t = 1
add!(ch, Dirichlet(:u, getfacetset(grid, "right"), f_dbc, 1))
close!(ch);

# Define cellvalues 
qr = QuadratureRule{RefQuadrilateral}(2)
cv = CellValues(qr, ip);

# ## FerriteAssembly setup 
# We first define the material models, and use the `ReducedStressState` to get a 
# plane strain response. 
elastic_material = ReducedStressState(PlaneStrain(), LinearElastic(;E=210e3, ν=0.3))

plastic_material = ReducedStressState(PlaneStrain(), Plastic(;
    elastic   = LinearElastic(E = 210e3, ν = 0.3),
    yield     = 100.0, 
    isotropic = Voce(;Hiso = 100e3, κ∞ = 1000.0),
    kinematic = ArmstrongFrederick(;Hkin = 0.0, β∞ = 1.0) # No kinematic hardening
    ));

# We need to create the domain buffers, where the difference from earlier tutorials 
# is that we have multiple domains. The domain specification should then be a `Dict`
# with one entry for each domain:
domains = Dict(
    "elastic"=>DomainSpec(dh, elastic_material, cv; set=getcellset(grid, "inclusion")),
    "plastic"=>DomainSpec(dh, plastic_material, cv; set=getcellset(grid, "matrix")) );

# Now, we can call `setup_domainbuffers`, which accepts the same keyword arguments 
# as `setup_domainbuffer`. Here, we accept the defaults. 
buffer = setup_domainbuffers(domains);

# ### Postprocessing setup 
# In this tutorial, we also demonstrate how the `QuadratureEvaluator` can be used 
# to obtain quadrature point data which can be used to visualize the results, in this 
# case the stresses. Specifically, we will use the evaluated data in combination with 
# `Ferrite`'s `L2Projector`.
# 
# First, we define a function to calculate the stresses for each material.
# Note that here we have to use some internals from `MechanicalMaterialModels.jl`,
# but this should be solved with 
# [MaterialModelsBase#12](https://github.com/KnutAM/MaterialModelsBase.jl/issues/12).

function calculate_stress(m::ReducedStressState, u, ∇u, qp_state)
    ϵ = MaterialModelsBase.expand_tensordim(m.stress_state, symmetric(∇u))
    σ = calculate_stress(m.material, ϵ, qp_state)
    return MaterialModelsBase.reduce_tensordim(m.stress_state, σ)
end
calculate_stress(m::LinearElastic, ϵ, qp_state) = m.C ⊡ ϵ
calculate_stress(m::Plastic, ϵ, qp_state) = calculate_stress(m.elastic, ϵ - qp_state.ϵp, qp_state);

# And then we create the QuadratureEvaluator including this function
qe = QuadratureEvaluator{SymmetricTensor{2,2,Float64,3}}(buffer, calculate_stress);

# Finally, we'll setup the L2Projector that we will use
proj = L2Projector(grid)
add!(proj, 1:getncells(grid), ip; qr_rhs = qr)
close!(proj);

# ## Solving the nonlinear problem via time-stepping
function solve_nonlinear_timehistory(buffer, dh, ch, l2_proj, qp_evaluator; time_history)
    maxiter = 100
    tolerance = 1e-6
    K = allocate_matrix(dh)
    r = zeros(ndofs(dh))
    a = zeros(ndofs(dh))
    ## Prepare postprocessing
    pvd = paraview_collection("multiple_materials")
    for (n, t) in enumerate(time_history)
        ## Update and apply the Dirichlet boundary conditions
        update!(ch, t)
        apply!(a, ch)
        for i in 1:maxiter
            ## Assemble the system 
            assembler = start_assemble(K, r)
            work!(assembler, buffer; a=a)
            ## Apply boundary conditions
            apply_zero!(K, r, ch)
            ## Check convergence
            norm(r) < tolerance && break
            i == maxiter && error("Did not converge")
            ## Solve the linear system and update the dof vector
            a .-= K \ r
            apply!(a, ch)
        end
        
        ## If converged, update the old state variables to the current. 
        update_states!(buffer)

        ## Postprocess
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
solve_nonlinear_timehistory(buffer, dh, ch, proj, qe; time_history=collect(range(0, 1, 20)));

#md # ## [Plain program](@id mixed_materials_plain_program)
#md #
#md # Here follows a version of the program without any comments.
#md # The file is also available here: [`mixed_materials.jl`](mixed_materials.jl).
#md #
#md # ```julia
#md # @__CODE__
#md # ```