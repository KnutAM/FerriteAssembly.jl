# Workers
The to perform work, such as assembling, the function [`work!`](@ref) should be called with a 
given worker and domain buffer. 

There are currently two categories of workers implemented:  which are used to calculate
- [Assemblers](@ref System-matrix-and-vector): Calculate system matrices and vectors
- [Integrators](@ref Integration): Integrate a function over the domain

```@docs
work!
```

## System matrix and vector
In order to assemble the system matrices and residual vectors, it is necessary to define the appropriate 
element routine and choose an assembler. The following assemblers can be used to assemble the system matrix and vector:
- [`Ferrite.start_assemble`](https://ferrite-fem.github.io/Ferrite.jl/stable/reference/assembly/#Ferrite.start_assemble) that returns
  - `Ferrite.AssemblerSparsityPattern`
  - `Ferrite.AssemblerSymmetricSparsityPattern`
- [`KeReAssembler`](@ref)
- [`ReAssembler`](@ref)

Where the ones available in `FerriteAssembly` have additional features, 
such as the possibility of applying constraints locally (see `Ferrite.apply_assemble!`) or calculate [scaling](@ref ElementResidualScaling) for the residual. 

### Element routines
One of the element methods should be overloaded for `material`. 
Note that `cellvalues` are already `reinit!`:ed for the current cell.
```@docs
FerriteAssembly.element_routine!(args...)
FerriteAssembly.element_residual!
```

### Assemblers
```@docs
KeReAssembler
ReAssembler
```

## Integration
In addition to assembling system matrices and vectors, much of the internal code can be reused 
to create quite efficient integration of values, given a solution vector (and potentially state variables).
The same domain buffer as for assembly can be used, and exactly as for assembly, we just call 
`work!(integrator, buffer)`. The following integrators are implemented
```@docs
SimpleIntegrator
Integrator
```
### Integration routines
The general integrator requires overloading the `integrate_cell!` function
```@docs
FerriteAssembly.integrate_cell!
```