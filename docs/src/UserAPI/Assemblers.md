# Assemblers
The assembly process runs by calling [`doassemble!`](@ref) with a given assembler. 
There are currently two categories of assemblers implemented, which are used to calculate
- [System matrices and vectors](@ref Assembling-system-matrices-and-residual-vectors)
- [Integrating a function](@ref Integration)

```@docs
doassemble!
```

## Assembling system matrices and residual vectors
In order to assemble the system matrices and residual vectors, it is necessary to define the appropriate 
element routine and choose an assembler. The following assemblers can be used to assemble the system matrix and vector:
- `Ferrite.AssemblerSparsityPattern` or `Ferrite.AssemblerSymmetricSparsityPattern` created with [`Ferrite.start_assemble`](https://ferrite-fem.github.io/Ferrite.jl/stable/reference/assembly/#Ferrite.start_assemble)
- [`KeReAssembler`](@ref)
- [`ReAssembler`](@ref)

Where the ones available in `FerriteAssembly` have additional features, 
such as the possibility of applying constraints locally (see `Ferrite.apply_assemble`) or calculate [scaling](@ref ElementResidualScaling) for the residual. 

### Element routines
One of the element methods should be overloaded for `material`. 
Note that `cellvalues` are already `reinit!`:ed for the current cell.
```@docs
FerriteAssembly.element_routine!
FerriteAssembly.element_residual!
```

### Assemblers
```@docs
KeReAssembler
ReAssembler
```

## Integration
In addition to assembling system matrices and vectors, much of the internal code can be reused 
to create quite efficient integration of values, given a solution vector (and potentially state variables)
The general workflow assumes that [`setup_assembly`](@ref) has already been called, and that a solution 
vector is available. Then, it is possible to call [`doassemble!`](@ref) with an integrator to 
obtain the integrated value. The following integrators are implemented
```@docs
SimpleIntegrator
Integrator
```