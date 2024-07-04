# # Multiple materials
# In this tutorial, we will see how we can assemble a domain and solve a problem 
# where we have multiple material behaviors. In this simple case, we will consider 
# an elastic inclusion, embedded in a plastically deforming matrix. 
# We will use models already implemented in the `ExampleElements` module, specifically 
# the `J2Plasticity` model (modified version of the model in the Ferrite plasticity 
# [example](https://ferrite-fem.github.io/Ferrite.jl/v0.3.14/examples/plasticity/), 
# adapted to follow the `MaterialModelsBase` interface) and the standard 
# `LinearElastic` model. 
# 
# The full script without intermediate comments is available at the 
# [bottom of this page](@ref mixed_materials_plain_program).
#
using Ferrite, FerriteAssembly, MaterialModelsBase
import FerriteAssembly.ExampleElements: J2Plasticity, ElasticPlaneStrain

# ## Setup Ferrite quantities
# We start by the grid with sets for the inclusion with 
# radius 0.8 and the surrounding matrix. 
function create_grid_with_inclusion()
    p1 = Vec((-1.0, -1.0))
    p2 = Vec(( 1.0,  1.0))
    grid = generate_grid(Quadrilateral, (20,20), p1, p2)
    addcellset!(grid, "inclusion", x -> norm(x) < 0.8)
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
add!(ch, Dirichlet(:u, getfacetset(grid,"left"), Returns(zero(Vec{2}))))
f_dbc(x,t) = Vec((0.05*t, 0.0))
add!(ch, Dirichlet(:u, getfacetset(grid, "right"), f_dbc))
close!(ch);

# Define cellvalues 
qr = QuadratureRule{RefQuadrilateral}(2)
cv = CellValues(qr, ip);

# ## FerriteAssembly setup 
# We first define materials for a plane strain state. The `J2Plasticity` model 
# is only implemented for 3d, so we use the `ReducedStressState` wrapper from 
# `MaterialModelsBase.jl` together with the `PlaneStrain` state. 
elastic_material = ElasticPlaneStrain(;E=210e3, ν=0.3)
plastic_material = ReducedStressState(
    PlaneStrain(), 
    J2Plasticity(;E=210e3, ν=0.3, σ0=100.0, H=10e3));

# We need to create the domain buffers, where the difference from earlier tutorials 
# is that we have multiple domains. The domain specification should then be a `Dict`
# with one entry for each domain:
domains = Dict(
    "elastic"=>DomainSpec(dh, elastic_material, cv; set=getcellset(grid, "inclusion")),
    "plastic"=>DomainSpec(dh, plastic_material, cv; set=getcellset(grid, "matrix")) );

# Now, we can call `setup_domainbuffers`, which accepts the same keyword arguments 
# as `setup_domainbuffer`. Here, we accept the defaults. 
buffer = setup_domainbuffers(domains);

# ## Solving the nonlinear problem via time-stepping
function solve_nonlinear_timehistory(buffer, dh, ch; time_history)
    maxiter = 10
    tolerance = 1e-6
    K = allocate_matrix(dh)
    r = zeros(ndofs(dh))
    a = zeros(ndofs(dh))
    ## Prepare postprocessing
    pvd = VTKFileCollection("multiple_materials.pvd", dh);
    addstep!(pvd, 0.0) do io
        write_solution(io, dh, a)
        Ferrite.write_cellset(io, dh.grid, "inclusion")
    end
    for t in time_history
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
            a .-= K\r
            apply!(a, ch)
        end
        ## Postprocess
        addstep!(pvd, 0.0) do io
            write_solution(io, dh, a)
            Ferrite.write_cellset(io, dh.grid, "inclusion")
        end
        ## If converged, update the old state variables to the current. 
        update_states!(buffer)
    end
    close(pvd)
    return nothing
end;
solve_nonlinear_timehistory(buffer, dh, ch; time_history=collect(range(0,1,20)));

#md # ## [Plain program](@id mixed_materials_plain_program)
#md #
#md # Here follows a version of the program without any comments.
#md # The file is also available here: [`mixed_materials.jl`](mixed_materials.jl).
#md #
#md # ```julia
#md # @__CODE__
#md # ```