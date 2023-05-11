# The full example can be downloaded [here](firstexample.jl).
# 
# First we create the dofhandler and cellvalues as in 
# [`Ferrite.jl`'s heat equation example](https://ferrite-fem.github.io/Ferrite.jl/stable/examples/heat_equation/)
using Ferrite, FerriteAssembly, BenchmarkTools
dh = DofHandler(generate_grid(Quadrilateral, (100, 100))); add!(dh, :u, 1); close!(dh)
cellvalues = CellScalarValues(QuadratureRule{2, RefCube}(2), Lagrange{2, RefCube, 1}());

# We start by defining the material 
# (that normally contains material parameters but are hard-coded in the example)
struct ThermalMaterial end;

# and then define our `element_routine!` for that material as
function FerriteAssembly.element_routine!(
        Ke::AbstractMatrix, re::AbstractVector, state, ae::AbstractVector, 
        material::ThermalMaterial, cellvalues::CellScalarValues, buffer::FerriteAssembly.CellBuffer
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
end;
# which is basically the same as in `Ferrite.jl`'s example. 

# We first call `setup_assembly` to setup all the variables required to assemble:
buffer, old_states, new_states = setup_assembly(dh, ThermalMaterial(), cellvalues);

# We can now create the global residual vectors and stiffness matrix
K = create_sparsity_pattern(dh)
r = zeros(ndofs(dh));

# And then assemble them.
assembler = start_assemble(K, r)
doassemble!(assembler, new_states, buffer);
K1 = deepcopy(K); #hide
# where we only need to give the `new_states`, since we do not rely on the old states.

# ### Threaded assembly
# To do the assembly in the example above threaded, we just tell `setup_assembly` that:
buffer2, _, _ = setup_assembly(dh, ThermalMaterial(), cellvalues; threading=true);
# This creates a default coloring of the grid, but custom coloring can also be given.

# We can call `doassemble!` as before
assembler = start_assemble(K, r)
doassemble!(assembler, new_states, buffer2);
K2 = deepcopy(K) #hide

# ### Automatic differentiation
# For elements for nonlinear coupled behavior, it can be time-consuming 
# to derive the analytical element stiffness. In these cases, it can be 
# advantageous to use automatic differentiation instead. This can be done 
# automatically by defining [`element_residual!`](@ref) instead of 
# [`element_routine!`](@ref) for the given material.
struct ThermalMaterialAD end

function FerriteAssembly.element_residual!(
        re::AbstractVector, state, ae::AbstractVector, 
        material::ThermalMaterialAD, cellvalues, buffer
        )
    n_basefuncs = getnbasefunctions(cellvalues)
    ## Loop over quadrature points
    for q_point in 1:getnquadpoints(cellvalues)
        ## Get the quadrature weight
        dΩ = getdetJdV(cellvalues, q_point)
        ∇u = function_gradient(cellvalues, q_point, ae)
        ## Loop over test shape functions
        for i in 1:n_basefuncs
            δN  = shape_value(cellvalues, q_point, i)
            ∇δN = shape_gradient(cellvalues, q_point, i)
            ## re = fint - fext
            re[i] += (∇δN ⋅ ∇u - δN) * dΩ
        end
    end
end;

buffer_ad, old_states_ad, new_states_ad = setup_assembly(dh, ThermalMaterialAD(), cellvalues);

# In this case we need the `ae` input and must therefore define `a`:
a = zeros(ndofs(dh))
assembler = start_assemble(K, r)
doassemble!(assembler, new_states_ad, buffer_ad; a=a);
K3 = deepcopy(K) #hide

# However, explicitly defining the element stiffness was a lot faster and has less allocations
@btime doassemble!($assembler, $new_states, $buffer; a=$a)
@btime doassemble!($assembler, $new_states_ad, $buffer_ad; a=$a)

# By using the special [`FerriteAssembly.AutoDiffCellBuffer`](@ref) that caches some variables for
# automatic differentiation, we can significantly improve performance. 
buffer_ad2, _, _ = setup_assembly(dh, ThermalMaterialAD(), cellvalues; autodiffbuffer=true)
@btime doassemble!($assembler, $new_states_ad, $buffer_ad2; a=$a)

assembler = start_assemble(K, r)                            #hide
doassemble!(assembler, new_states_ad, buffer_ad2; a=a);     #hide
K4 = deepcopy(K)                                            #hide

# Test that all methods give the same stiffness #src
using Test                                      #hide
@test K1 ≈ K2 ≈ K3 ≈ K4                         #hide