# The full example can be downloaded [here](firstexample.jl).
# In order to use FerriteAssembly for assembly, we need to 
# 1) Define a custom `material` that determines the physics 
# 2) Overload `element_routine!` or `element_residual!` for that material.
# 
# We can then set up Ferrite in the standard way, by defining a 
# Grid, DofHandler, and CellValues. We then pass these to the 
# setup_domainbuffer function and create an assembler (worker),
# allowing us to perform the work with the `work!` function.
# 
# ### Standard Ferrite setup
# First we create the dofhandler and cellvalues as in 
# [`Ferrite.jl`'s heat equation example](https://ferrite-fem.github.io/Ferrite.jl/stable/examples/heat_equation/)
using Ferrite, FerriteAssembly, BenchmarkTools
dh = DofHandler(generate_grid(Quadrilateral, (100, 100))); add!(dh, :u, 1); close!(dh)
cellvalues = CellScalarValues(QuadratureRule{2, RefCube}(2), Lagrange{2, RefCube, 1}());
K = create_sparsity_pattern(dh)
r = zeros(ndofs(dh));

# ### Define the physics
# We start by defining the material and create an instance of it
struct ThermalMaterial 
    k::Float64 # Thermal conductivity
    f::Float64 # Volumetric heat source
end
material = ThermalMaterial(1.0, 1.0)

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

# ### Setting up the domain.
# We first start by defining a domain and passing 
# that to the `setup_domainbuffer` function. 
grid_domain = DomainSpec(dh, material, cellvalues)
buffer = setup_domainbuffer(grid_domain);

# ### Performing the assembly
# The `worker` in this case is the standard Ferrite assembler:
assembler = start_assemble(K, r);

# Given this worker, we can do the work to assemble `K` and `r`
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

# ### Automatic Differentiation
# Deriving the analytical element stiffness is often time-consuming,
# and can be avoided by using Automatic Differentiation (AD). 
# FerriteAssembly supports this: Just define [`element_residual!`](@ref) 
# instead of [`element_routine!`](@ref) for your material.
struct ThermalMaterialAD
    k::Float64 # Thermal conductivity
    f::Float64 # Volumetric heat source
end
material_ad = ThermalMaterialAD(1.0, 1.0)

function FerriteAssembly.element_residual!(re, state, ae, 
        material::ThermalMaterialAD, cellvalues, cellbuffer
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

# We use the same setup as before, but with `material_ad`
grid_domain_ad = DomainSpec(dh, material_ad, cellvalues)
buffer_ad = setup_domainbuffer(grid_domain_ad);

# In this case we need the `ae` input and must therefore define `a`:
a = zeros(ndofs(dh))
assembler = start_assemble(K, r)
work!(assembler, buffer_ad; a=a);
K3 = deepcopy(K); #hide

# ### AD performance 
# Explicitly defining the element stiffness was a lot faster and has less allocations:
@btime work!($assembler, $buffer; a=$a)
@btime work!($assembler, $buffer_ad; a=$a)

# FerriteAssembly comes with a special `cellbuffer` for speeding up 
# automatic differentiation, we can significantly improve the performance.
buffer_ad2 = setup_domainbuffer(grid_domain_ad; autodiffbuffer=true)
@btime work!($assembler, $buffer_ad2; a=$a)
                                                            #hide
assembler = start_assemble(K, r)                            #hide
work!(assembler, buffer_ad2; a=a);                          #hide
K4 = deepcopy(K);                                           #hide
                                                            #hide
# Test that all methods give the same stiffness #src
using Test                                      #hide
@test K1 ≈ K2 ≈ K3 ≈ K4                         #hide
nothing;                                        #hide