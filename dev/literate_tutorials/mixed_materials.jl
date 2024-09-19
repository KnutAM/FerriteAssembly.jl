# # Multiple materials
# In this tutorial, we will see how we can assemble a domain and solve a problem 
# where we have multiple material behaviors. In this simple case, we will consider 
# an elastic inclusion, embedded in a plastically deforming matrix. 
# 
# ![results](mixed_materials.png)
# 
# **Figure 1:** Results showing horizontal stresses ($\sigma_{11}$ [MPa]) on the deformed mesh.
# 
# For this example, we'll use material models defined in the 
# [`MechanicalMaterialModels.jl`](https://github.com/KnutAM/MechanicalMaterialModels.jl)
# package, which defines models according to the `MaterialModelsBase` interface.
# 
# The full script without intermediate comments is available at the 
# [bottom of this page](@ref mixed_materials_plain_program).
#
using Ferrite, FerriteAssembly, FerriteMeshParser
using MaterialModelsBase, MechanicalMaterialModels, WriteVTK

# ## Setup Ferrite quantities
# We start by the loading a grid containing a central inclusion,
grid = get_ferrite_grid(joinpath(@__DIR__, "square_with_inclusion.inp"));

# Define interpolation
ip = Lagrange{RefTriangle, 2}()^2;

# Followed by the dof handler 
dh = DofHandler(grid)
add!(dh, :u, ip)
close!(dh);

# Dirichlet boundary conditions
ch = ConstraintHandler(dh)
add!(ch, Dirichlet(:u, getfacetset(grid,"left"), Returns(0.0), 1))
add!(ch, Dirichlet(:u, getfacetset(grid,"bottom"), Returns(0.0), 2))
close!(ch);

# and a Neumann boundary condition.
lh = LoadHandler(dh)
add!(lh, Neumann(:u, 3, getfacetset(grid, "right"), (x, t, n) -> 1e3 * t * n));

# Finally, we define the cellvalues 
qr = QuadratureRule{RefTriangle}(4)
ip_geo = Lagrange{RefTriangle, 2}()
cv = CellValues(qr, ip, ip_geo);

# ## FerriteAssembly setup 
# We first define the material models, and use the `ReducedStressState` to get a 
# plane stress response. 
elastic_material = ReducedStressState(PlaneStress(), LinearElastic(;E=210e3, ν=0.3))

plastic_material = ReducedStressState(PlaneStress(), Plastic(;
    elastic   = LinearElastic(E = 210e3, ν = 0.3),
    yield     = 100.0, 
    isotropic = Voce(;Hiso = 50e3, κ∞ = 1e9),         # Linear isotropic hardening
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
# In this tutorial, we also demonstrate how the `QuadPointEvaluator` can be used 
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

# And then we create the QuadPointEvaluator including this function
qe = QuadPointEvaluator{SymmetricTensor{2,2,Float64,3}}(buffer, calculate_stress);

# Finally, we'll setup the L2Projector that we will use
proj = L2Projector(grid)
add!(proj, 1:getncells(grid), Lagrange{RefTriangle, 1}(); qr_rhs = qr)
close!(proj);

# ## Solving the nonlinear problem via time-stepping
function solve_nonlinear_timehistory(buffer, dh, ch, lh, l2_proj, qp_evaluator; time_history)
    maxiter = 100
    tolerance = 1e-6
    K = allocate_matrix(dh)
    r = zeros(ndofs(dh))
    fext = zeros(ndofs(dh))
    a = zeros(ndofs(dh))
    ## Prepare postprocessing
    pvd = paraview_collection("multiple_materials")
    for (n, t) in enumerate(time_history)
        ## Update and apply the Dirichlet boundary conditions
        update!(ch, t)
        apply!(a, ch)
        apply!(fext, lh, t)
        for i in 1:maxiter
            ## Assemble the system
            assembler = start_assemble(K, r)
            work!(assembler, buffer; a=a)
            ## Apply boundary conditions
            r .-= fext
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
    vtk_save(pvd)
    return nothing
end;
solve_nonlinear_timehistory(buffer, dh, ch, lh, proj, qe; time_history=collect(range(0, 1, 20)));

#md # ## [Plain program](@id mixed_materials_plain_program)
#md #
#md # Here follows a version of the program without any comments.
#md # The file is also available here: [`mixed_materials.jl`](mixed_materials.jl).
#md #
#md # ```julia
#md # @__CODE__
#md # ```
