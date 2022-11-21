# The full example can be downloaded [here](firstexample.jl).
# 
# First we create the dofhandler and cellvalues as in 
# [`Ferrite.jl`'s heat equation example](https://ferrite-fem.github.io/Ferrite.jl/stable/examples/heat_equation/)
using Ferrite, FerriteAssembly, BenchmarkTools
dh = DofHandler(generate_grid(Quadrilateral, (100, 100))); push!(dh, :u, 1); close!(dh)
cellvalues = CellScalarValues(QuadratureRule{2, RefCube}(2), Lagrange{2, RefCube, 1}())

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
            δN  = shape_value(cellvalues, q_point, i)
            ∇δN = shape_gradient(cellvalues, q_point, i)
            ## Add body load contribution to re
            re[i] += -δN * dΩ
            ## Loop over trial shape functions
            for j in 1:n_basefuncs
                ∇N = shape_gradient(cellvalues, q_point, j)
                ## Add contribution to Ke
                Ke[i, j] += (∇δN ⋅ ∇N) * dΩ
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

# The next step is gathering all variables that are shared for all elements,
# but may be modified for or by each element in the `cellbuffer`:
cellbuffer = CellBuffer(dh, cellvalues, ThermalMaterial())
# Note that `cellvalues` can be a `Tuple` or `NamedTuple`, 
# this is useful for coupled problems with multiple fields. 

# We then define our `assembler` and do the assembly that will modify `K` and `r`
assembler = start_assemble(K,r)
doassemble!(assembler, cellbuffer, states, dh)

# ### Threaded assembly
# To do the assembly in the example above threaded, 
# we need to color the grid to avoid race conditions.
# This can be done with `Ferrite.jl`'s `create_coloring` function:
colors = create_coloring(dh.grid)

# We must also create threaded versions of the `cellbuffer` and `assembler`,
cellbuffers = create_threaded_CellBuffers(CellBuffer(dh, cellvalues, ThermalMaterial()))
assemblers = create_threaded_assemblers(K, r)

# And then we can call `doassemble!` as
doassemble!(assemblers, cellbuffers, states, dh, colors)

# ### Automatic differentiation
# For elements for nonlinear coupled behavior, it can be time-consuming 
# to derive the analytical element stiffness. In these cases, it can be 
# advantageous to use automatic differentiation instead. This can be done 
# automatically by defining [`element_residual!`](@ref) instead of 
# [`element_routine!`](@ref) for the given material.
struct ThermalMaterialAD end 

function FerriteAssembly.element_residual!(
    re::AbstractVector, state,
    ae::AbstractVector, material::ThermalMaterialAD, cellvalues, 
    dh_fh::Union{DofHandler,FieldHandler}, Δt, buffer
    )
    n_basefuncs = getnbasefunctions(cellvalues)
    # Loop over quadrature points
    for q_point in 1:getnquadpoints(cellvalues)
        # Get the quadrature weight
        dΩ = getdetJdV(cellvalues, q_point)
        ∇u = function_gradient(cellvalues, q_point, ae)
        # Loop over test shape functions
        for i in 1:n_basefuncs
            δN  = shape_value(cellvalues, q_point, i)
            ∇δN = shape_gradient(cellvalues, q_point, i)
            # re = fint - fext
            re[i] += (∇δN ⋅ ∇u - δN) * dΩ
        end
    end
end

# In this case we need the `ae` input and must therefore define 
a = zeros(ndofs(dh))
cellbuffer2 = CellBuffer(dh, cellvalues, ThermalMaterialAD())
assembler = start_assemble(K,r)
doassemble!(assembler, cellbuffer2, states, dh, a)

# However, explicitly defining the element stiffness was a lot faster and has less allocations
@btime doassemble!(assembler, $cellbuffer, $states, $dh) setup=(assembler=start_assemble(K,r))
@btime doassemble!(assembler, $cellbuffer2, $states, $dh) setup=(assembler=start_assemble(K,r))

# But by defining a special cellbuffer that caches some variables for use with automatic differentiation,
# we can get closer to the performance of explicitly defining the element stiffness
cellbuffer3 = FerriteAssembly.AutoDiffCellBuffer(states, dh, cellvalues, ThermalMaterialAD())
@btime doassemble!(assembler, $cellbuffer3, $states, $dh) setup=(assembler=start_assemble(K,r))