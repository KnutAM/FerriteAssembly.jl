# Assemblers
```@docs
doassemble!
```

## Assembling system matrices and residual vectors


### Element routines
One of the element methods should be overloaded for `material`. 
Note that `cellvalues` are already `reinit!`:ed for the current cell.
```@docs
FerriteAssembly.element_routine!
FerriteAssembly.element_residual!
```

### Assemblers
`Ferrite.AssemblerSparsityPattern` and `Ferrite.AssemblerSymmetricSparsityPattern` 
are supported. These are created as described in `Ferrite.jl`'s documentation, and
when passed to `doassemble!`, they assemble `K` and `r` as in `Ferrite.jl`. In addition,
the following assemblers are provided by `FerriteAssemble.jl`,
```@docs
KeReAssembler
ReAssembler
```

## Integration
