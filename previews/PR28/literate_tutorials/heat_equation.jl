# # Heat Equation
# This tutorial solves the stationary heat flow example, 
# equivalent to the first example in `Ferrite.jl`. 
# 
# The full script without intermediate comments is available at the 
# [bottom of this page](@ref heat_equation_plain_program).
# 
# ## Setup
# First we create the dofhandler, vectors and matrices, and cellvalues as in 
# [`Ferrite.jl`'s heat equation example](https://ferrite-fem.github.io/Ferrite.jl/stable/examples/heat_equation/)
using Ferrite, FerriteAssembly
grid = generate_grid(Quadrilateral, (20, 20))
ip = Lagrange{RefQuadrilateral,1}()
dh = DofHandler(grid); add!(dh, :u, ip); close!(dh)
cellvalues = CellValues(QuadratureRule{RefQuadrilateral}(2), ip);
K = allocate_matrix(dh)
r = zeros(ndofs(dh));

# ## Define the physics
# We start by defining the material and create an instance of it
struct ThermalMaterial 
    k::Float64 # Thermal conductivity
    f::Float64 # Volumetric heat source
end
material = ThermalMaterial(1.0, 1.0);

# and then define our `element_routine!` for that material as
function FerriteAssembly.element_routine!(Ke, re, state, ae, 
        material::ThermalMaterial, cellvalues, cellbuffer
        )
    n_basefuncs = getnbasefunctions(cellvalues)
    ## Loop over quadrature points
    for q_point in 1:getnquadpoints(cellvalues)
        dΩ = getdetJdV(cellvalues, q_point)
        for i in 1:n_basefuncs
            δN  = shape_value(cellvalues, q_point, i)
            ∇δN = shape_gradient(cellvalues, q_point, i)
            ## Add body load contribution to re
            re[i] += -material.f*δN * dΩ
            ## Loop over trial shape functions
            for j in 1:n_basefuncs
                ∇N = shape_gradient(cellvalues, q_point, j)
                ## Add contribution to Ke
                Ke[i, j] += material.k*(∇δN ⋅ ∇N) * dΩ
            end
        end
    end
end;
# which is basically the same as in `Ferrite.jl`'s example. 

# ## Assemble 
# We first start by defining a domain and passing 
# that to the `setup_domainbuffer` function. 
grid_domain = DomainSpec(dh, material, cellvalues)
buffer = setup_domainbuffer(grid_domain);

# The `worker` in this case is the standard Ferrite assembler:
assembler = start_assemble(K, r);

# Given this worker, we can do the work to assemble `K` and `r`
work!(assembler, buffer);

# ## Solve the problem. 
# To actually solve the problem, we also need Dirichlet boundary conditions.
ch = ConstraintHandler(dh)
facetset = union((getfacetset(grid,k) for k in ("left", "right", "bottom", "top"))...)
add!(ch, Dirichlet(:u, facetset, Returns(0.0)))
close!(ch);
apply_zero!(K, r, ch)
# where we use `apply_zero!` since we assembled assuming a zero temperature. 
# We could have used `apply!(K,f,ch)`, but for non-zero dirichlet conditions, 
# this relies on the correct sign for external loads, and we have r=-f.

# Finally, we can solve the problem and save the results 
a = -K\r
VTKFile("heat_equation", grid) do vtk
    write_solution(vtk, dh, a)
end;

#md # ## [Plain program](@id heat_equation_plain_program)
#md #
#md # Here follows a version of the program without any comments.
#md # The file is also available here: [`heat_equation.jl`](heat_equation.jl).
#md #
#md # ```julia
#md # @__CODE__
#md # ```
