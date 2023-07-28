# The full example can be downloaded [here](firstexample.jl).
# This example covers the following
# 1) How to assemble the heat equation in the standard way 
# 2) How to do the assembly threaded 
# 3) How to use automatic differentiation
# 
# First we create the dofhandler and cellvalues as in 
# [`Ferrite.jl`'s heat equation example](https://ferrite-fem.github.io/Ferrite.jl/stable/examples/heat_equation/)
using Ferrite, FerriteAssembly, BenchmarkTools
dh = DofHandler(generate_grid(Quadrilateral, (100, 100))); add!(dh, :u, 1); close!(dh)
cellvalues = CellScalarValues(QuadratureRule{2, RefCube}(2), Lagrange{2, RefCube, 1}());

# We start by defining the material 
struct ThermalMaterial 
    k::Float64 # Thermal conductivity
    f::Float64 # Volumetric heat source
end

# and then define our `element_routine!` for that material as
function FerriteAssembly.element_routine!(Ke, re, state, ae, 
        material::ThermalMaterial, cellvalues, buffer
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

# We first call `setup_assembly` to setup all the variables required to assemble:
material = ThermalMaterial(1.0, 1.0)
grid_domain = GridDomain(dh, material, cellvalues)
buffer = setup_domainbuffer(grid_domain);

# We can now create the global residual vectors and stiffness matrix
K = create_sparsity_pattern(dh)
r = zeros(ndofs(dh));

# And then assemble them.
assembler = start_assemble(K, r)
work!(assembler, buffer);
K1 = deepcopy(K); #hide

# ### Threaded assembly
# To do the assembly in the example above threaded, we just tell `setup_domainbuffer` that:
threaded_buffer = setup_domainbuffer(grid_domain; threading=true);
# This creates a default coloring of the grid, but custom coloring can also be given.

# We can call `work!` as before
assembler = start_assemble(K, r)
work!(assembler, threaded_buffer);
K2 = deepcopy(K); #hide

# ### Automatic differentiation
# For elements for nonlinear coupled behavior, it can be time-consuming 
# to derive the analytical element stiffness. In these cases, it can be 
# advantageous to use automatic differentiation instead. This can be done 
# automatically by defining [`element_residual!`](@ref) instead of 
# [`element_routine!`](@ref) for the given material.
struct ThermalMaterialAD
    k::Float64 # Thermal conductivity
    f::Float64 # Volumetric heat source
end

function FerriteAssembly.element_residual!(re, state, ae, 
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
            re[i] += (material.k*(∇δN ⋅ ∇u) - material.f*δN) * dΩ
        end
    end
end;

# We then create an instance of this material, and setup the assembly,
material_ad = ThermalMaterialAD(1.0, 1.0)
grid_domain_ad = GridDomain(dh, material_ad, cellvalues)
buffer_ad = setup_domainbuffer(grid_domain_ad);

# In this case we need the `ae` input and must therefore define `a`:
a = zeros(ndofs(dh))
assembler = start_assemble(K, r)
work!(assembler, buffer_ad; anew=a);
K3 = deepcopy(K); #hide

# However, explicitly defining the element stiffness was a lot faster and has less allocations
@btime work!($assembler, $buffer; anew=$a)
@btime work!($assembler, $buffer_ad; anew=$a)

# By using the special [`FerriteAssembly.AutoDiffCellBuffer`](@ref) that caches some variables for
# automatic differentiation, we can significantly improve the performance.
buffer_ad2 = setup_domainbuffer(grid_domain_ad; autodiffbuffer=true)
@btime work!($assembler, $buffer_ad2; anew=$a)
                                                            #hide
assembler = start_assemble(K, r)                            #hide
work!(assembler, buffer_ad2; anew=a);                       #hide
K4 = deepcopy(K);                                           #hide
                                                            #hide
# Test that all methods give the same stiffness #src
using Test                                      #hide
@test K1 ≈ K2 ≈ K3 ≈ K4                         #hide
nothing;                                        #hide