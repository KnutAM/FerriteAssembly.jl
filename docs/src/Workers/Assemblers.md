# Assemblers
Assemblers are used to assemble the system matrices and residual vectors. To calculate the contribution to these matrices, one or more of the following functions should be overloaded for the specific `material`:
* [`element_routine!`](@ref FerriteAssembly.element_routine!)
* [`element_residual!`](@ref FerriteAssembly.element_residual!)
* [`face_residual!`](@ref FerriteAssembly.face_residual!)
* [`face_routine!`](@ref FerriteAssembly.face_routine!)


## Available assemblers 
The following assemblers can be used to assemble the system matrix and vector:
- [`Ferrite.start_assemble`](@ref Ferrite-assemblers) (`AssemblerSparsityPattern` or `AssemblerSymmetricSparsityPattern`)
- [`KeReAssembler`](@ref)
- [`ReAssembler`](@ref)

### Ferrite assemblers 
Ferrite assemblers are supported and are used to assemble both a matrix and a vector. They are created by calling Ferrite's `start_assemble` function. They will primarily call `element_routine!` or `face_routine!`. 
However, if not defined for the given material, automatic differentiation will be used to get the matrix contribution by calling `element_residual!` or `face_residual!` instead.

### ReAssembler
The `ReAssembler` assembles only the residual vector, and will thus call `element_residual!` or `face_residual!`.
```@docs 
ReAssembler
```

### KeReAssembler
The `KeReAssembler` assembles both the residual vector and the system matrix, similar to Ferrite's assemblers. It will also primarily call `element_routine!` or `face_routine!` and fall back to automatic differentiation of `element_residual!` and `face_residual!` if the former methods are not defined. 
```@docs
KeReAssembler
```

## Functions to overload
```@docs
FerriteAssembly.element_routine!
FerriteAssembly.element_residual!
FerriteAssembly.face_routine!
FerriteAssembly.face_residual!
```

## Residual scaling
There are many options for how to scale the residual in finite element simulations.
This package does not intend to implement many different options, but does give the 
user the option to calculate scaling contributions from each cell, which may be useful.
By defining a `scaling` that is passed to [`KeReAssembler`](@ref) or [`ReAssembler`](@ref), it can be updated based on the output from each cell. 

One type of scaling, [`ElementResidualScaling`](@ref), is included as described below. Its code can be used as a template for how to include custom scaling that works on the element level.

### ElementResidualScaling
```@docs
ElementResidualScaling
```
