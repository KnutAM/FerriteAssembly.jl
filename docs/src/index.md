```@meta
CurrentModule = FerriteAssembly
```

# FerriteAssembly

The goal of [FerriteAssembly](https://github.com/KnutAM/FerriteAssembly.jl) 
is to provide a simple structure to perform assembly in 
[Ferrite.jl](https://github.com/Ferrite-FEM/Ferrite.jl/).

Both the `DofHandler` and `MixedDofHandler` are supported.
Currently only sequential assembly, but this will be extended 
to threaded assembly.

The package work by exporting the `doassemble!` function, and requires the 
user to define either `element_routine!` (calculate both `Ke` and `re`),
or just `element_residual!` (calculate only `re`). 
In the latter case, `Ke` is calculated by 
[ForwardDiff.jl](https://github.com/JuliaDiff/ForwardDiff.jl)

An advanced option to scale the unknowns, residual, and jacobian exists 
(currently undocumented)

```@contents
```

## A first example
We take a slightly modified element routine from `Ferrite.jl`'s heat equation 
[example](https://ferrite-fem.github.io/Ferrite.jl/stable/examples/heat_equation/).
We assume that we already have defined `dh::DofHandler` and our 
`cellvalues::CellScalarValues` according to that example. We then define
```julia
struct ThermalMaterial end
material = ThermalMaterial()
```
and our `element_routine!` for that material as 
```julia
function FerriteAssembly.element_routine!(
    Ke::AbstractMatrix, re::AbstractVector, 
    ae_new::AbstractVector, ae_old::AbstractVector,
    state::AbstractVector, material::ThermalMaterial, 
    cellvalues::CellScalarValues, 
    dh_fh::Union{DofHandler,FieldHandler}, Δt, materialcache
    )
    n_basefuncs = getnbasefunctions(cellvalues)
    # Loop over quadrature points
    for q_point in 1:getnquadpoints(cellvalues)
        dΩ = getdetJdV(cellvalues, q_point)
        for i in 1:n_basefuncs
            δu  = shape_value(cellvalues, q_point, i)
            ∇δu = shape_gradient(cellvalues, q_point, i)
            # Add body load contribution to re
            re[i] += -δu * dΩ
            # Loop over trial shape functions
            for j in 1:n_basefuncs
                ∇u = shape_gradient(cellvalues, q_point, j)
                # Add contribution to Ke
                Ke[i, j] += (∇δu ⋅ ∇u) * dΩ
            end
        end
    end
end
```
We can now create the global vectors and matrices, as well as the time step
```julia
K = create_sparsity_pattern(dh)
a=zeros(ndofs(dh)); aold=copy(a);
r = zeros(ndofs(dh))
Δt=1.0
```
Before we can assemble, we need our `cellcache` (to avoid unecessary allocations 
inside the assembly loop). But this is created automatically from the DofHandler:
```julia
cache = CellCache(dh)
```

And then we can assemble our system by calling 
```julia
doassemble!(K, r, a, aold, states, dh, cv, material, Δt, cache)
```

## Detailed API description
One of the element methods should be overloaded for a given combination of `cellvalues`
and `material`.
```@docs
FerriteAssembly.element_routine!
FerriteAssembly.element_residual!
```

The allocations are reduced by saving all variables to a cellcache, 
which is automatically created by giving it the dof handler:
```@docs
CellCache
```

Once everything is set up, one can call the function which will actually 
do the assembly:
```@docs
doassemble!
```
