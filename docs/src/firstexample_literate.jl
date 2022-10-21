# ## A minimal example
# The full example can be downloaded [here](firstexample.jl).
# 
# First we create the dofhandler and cellvalues as in 
# [`Ferrite.jl`'s heat equation example](https://ferrite-fem.github.io/Ferrite.jl/stable/examples/heat_equation/)
using Ferrite, FerriteAssembly
dh = DofHandler(generate_grid(Quadrilateral, (20, 20))); push!(dh, :u, 1); close!(dh);
cellvalues = CellScalarValues(QuadratureRule{dim, RefCube}(2), Lagrange{dim, RefCube, 1}());

# We start by defining the material 
# (that normally contains material parameters but are hard-coded in the example)
struct ThermalMaterial end

# and then define our `element_routine!` for that material as
function FerriteAssembly.element_routine!(
    Ke::AbstractMatrix, re::AbstractVector, state, 
    ae::AbstractVector, material::ThermalMaterial, cellvalues::CellScalarValues, 
    dh_fh, Δt, buffer::CellBuffer
    )
    n_basefuncs = getnbasefunctions(cellvalues)
    ## Loop over quadrature points
    for q_point in 1:getnquadpoints(cellvalues)
        dΩ = getdetJdV(cellvalues, q_point)
        for i in 1:n_basefuncs
            δu  = shape_value(cellvalues, q_point, i)
            ∇δu = shape_gradient(cellvalues, q_point, i)
            ## Add body load contribution to re
            re[i] += -δu * dΩ
            ## Loop over trial shape functions
            for j in 1:n_basefuncs
                ∇u = shape_gradient(cellvalues, q_point, j)
                ## Add contribution to Ke
                Ke[i, j] += (∇δu ⋅ ∇u) * dΩ
            end
        end
    end
end
# which is basically the same as in `Ferrite.jl`'s example. 

# We can now create the global vectors and matrices.
K = create_sparsity_pattern(dh)
r = zeros(ndofs(dh))

# In general, each integration point (or cell if desired by the user) can have a `state`.
# In this case, we have no state variables, but the interface still requires them.
# We can then create states of `nothing` by 
states = create_states(dh)
# See [`create_states`](@ref) for more detailed options of creating actual state variables. 

# The next step is gathering all variables that are modified for each cell, 
# but don't belong to each cell, in the `cellbuffer`:
cellbuffer = CellBuffer(dh, cellvalues, ThermalMaterial())
# Note that `cellvalues` can be a `Tuple` or `NamedTuple`, this is useful for coupled 
# problems with multiple fields. 

# We then define our `assembler` and do the assembly
assembler = start_assemble(K,r)
doassemble!(assembler, cellbuffer, states, dh, a, aold, Δt)

# ## Threaded assembly
# To do the assembly in the example above threaded, 
# we need to color the grid to avoid race conditions.
# This can be done with `Ferrite.jl`'s `create_coloring` function:
colors = create_coloring(dh.grid)

# We must also create threaded versions of the `cellbuffer` and `assembler`,
cellbuffers = create_threaded_CellBuffers(CellBuffer(dh, cellvalues, ThermalMaterial()))
assemblers = create_threaded_assemblers(K, r)

# And then we can call `doassemble!` as
doassemble!(assemblers, cellbuffers, states, colors, dh)
