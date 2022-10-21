```@meta
CurrentModule = FerriteAssembly
```

# FerriteAssembly

The goal of [FerriteAssembly](https://github.com/KnutAM/FerriteAssembly.jl) 
is to provide a simple structure to perform assembly in 
[Ferrite.jl](https://github.com/Ferrite-FEM/Ferrite.jl/).

Sequential and threaded assembly of both the `DofHandler` and `MixedDofHandler` are supported.

The package work by exporting the `doassemble!` function, and requires the 
user to define either `element_routine!` (calculate both `Ke` and `re`),
or just `element_residual!` (calculate only `re`). 
In the latter case, `Ke` is calculated by 
[ForwardDiff.jl](https://github.com/JuliaDiff/ForwardDiff.jl)

An advanced option to scale the unknowns, residual, and jacobian exists 
(currently undocumented)

```@contents
```

## A minimal example
First we create the dofhandler and cellvalues as in 
[`Ferrite.jl`'s heat equation example](https://ferrite-fem.github.io/Ferrite.jl/stable/examples/heat_equation/)
```julia
dh = DofHandler(generate_grid(Quadrilateral, (20, 20))); push!(dh, :u, 1); close!(dh);
cellvalues = CellScalarValues(QuadratureRule{dim, RefCube}(2), Lagrange{dim, RefCube, 1}());
```
We start by defining the material
(that normally contains material parameters
but those are hard-coded in the example)
```julia
struct ThermalMaterial end
```
Then, we define our `element_routine!` for that material as
```julia
function FerriteAssembly.element_routine!(
    Ke::AbstractMatrix, re::AbstractVector, state, 
    ae::AbstractVector, material::ThermalMaterial, cellvalues::CellScalarValues, 
    dh_fh, Δt, buffer::CellBuffer
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
which is basically the same as in `Ferrite.jl`'s example. 

We can now create the global vectors and matrices.
```julia
K = create_sparsity_pattern(dh)
r = zeros(ndofs(dh))
```
In general, each integration point (or cell if desired by the user) can have a `state`.
In this case, we have no state variables, but the interface still requires them.
We can then create states of `nothing` by 
```julia
states = create_states(dh)
```
See [`create_states`](@ref) for more detailed options of creating actual state variables. 

The next step is gathering all variables that are modified for each cell, 
but don't belong to each cell, in the `cellbuffer`:
```julia
cellbuffer = CellBuffer(dh, cellvalues, ThermalMaterial())
```
Note that `cellvalues` can be a `Tuple` or `NamedTuple`, this is useful for coupled 
problems with multiple fields. 
We then define our `assembler` and do the assembly
```julia
assembler = start_assemble(K,r)
doassemble!(assembler, cellbuffer, states, dh, a, aold, Δt)
```

## Threaded assembly
To do the assembly in the example above threaded, 
we need to color the grid to avoid race conditions.
This can be done with `Ferrite.jl`'s `create_coloring` function:
```julia
colors = create_coloring(dh.grid)
```
We must also create threaded versions of the `cellbuffer` and `assembler`,
```julia
cellbuffers = create_threaded_CellBuffers(CellBuffer(dh, cellvalues, ThermalMaterial()))
assemblers = create_threaded_assemblers(K, r)
```

And then we can call `doassemble!` as
```julia
doassemble!(assemblers, cellbuffers, states, colors, dh)
```

## Detailed API description
One of the element methods should be overloaded for a given combination of `cellvalues`
and `material`.
```@docs
FerriteAssembly.element_routine!
FerriteAssembly.element_residual!
```

Variables that are used and modified for each cell of a certain type, 
but that don't belong to a specific cell, are collected in a `CellBuffer`.
```@docs
CellBuffer
```

For parallel assembly, we need a vector of `CellBuffer`s: 
One `CellBuffer` for each thread.
For the `MixedDofHandler`, we first loop over the type of cells,
so we need a tuple that contains a vector of `CellBuffer`s. 
Construction of this via `deepcopy` is implemented as 
```@docs
create_threaded_CellBuffers
```

Similarily, we need a vector of assemblers that is convieniently 
created by calling 
```@docs
create_threaded_assemblers
```

The initial state variables may vary depending on the position in the grid.
Furthermore, the datastructure depends on the type of dof handler, so
a convenience function exists that creates the correct variable
```@docs
create_states
```

Once everything is set up, one can call the function which will actually 
do the assembly:
```@docs
doassemble!
```
