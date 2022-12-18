```@meta
CurrentModule = FerriteAssembly
```

# API
Methods that should be overloaded and exported functions are described on this page.
It may also be useful to checkout the [Data structures](@ref).

## Element routines
One of the element methods should be overloaded for a given combination of `cellvalues`
and `material`. 
Note that the `cellvalues` are already `reinit!`:ed when passed to the element routines. 
```@docs
FerriteAssembly.element_routine!
FerriteAssembly.element_residual!
```

## CellBuffer
Variables that are used and modified for each cell of a certain type, 
but that don't belong to a specific cell, are collected in a `CellBuffer`.
```@docs
CellBuffer
setup_cellbuffer
AutoDiffCellBuffer
setup_ad_cellbuffer
getCellBuffer
```

### Access functions
The following access functions can be used to extract information from the 
`CellBuffer`
```@docs
Ferrite.getcoordinates
FerriteAssembly.get_aeold
FerriteAssembly.get_load
FerriteAssembly.get_cache
```

## State variables
The initial state variables may vary depending on the position in the grid.
Furthermore, the datastructure depends on the type of dof handler, so
a convenience function exists that creates the correct variable
```@docs
create_states
```

## `doassemble!`
Once everything is set up, one can call the function which will actually 
do the assembly:
```@docs
doassemble!
```

## Threaded assembly
For parallel assembly, we need a vector of `CellBuffer`s and assemblers, 
one for each thread. 

For the `MixedDofHandler`, the outer loop is over the type of cells,
so we need a tuple that contains vectors of `CellBuffer`s. 
Construction of this via `deepcopy` is implemented as
```@docs
create_threaded_CellBuffers
```

A vector of assemblers that is convieniently created by calling 
```@docs
create_threaded_assemblers
```

## Residual scaling
There are many options for how to scale the residual in finite element simulations.
This package does not intend to implement many different options, but does give the 
user the option to calculate scaling contributions from each cell, which may be useful.
By defining a `scaling` that is passed to [`doassemble!`](@ref), this can be updated in each cell. 
One type of scaling, [`ElementResidualScaling`](@ref), is included as described below.
Its code can be used as a template 
for how to include custom scaling that works on the element level. 

In general, a scaling must support the [`update_scaling!`](@ref) function.
For consistency it is also nice, but not required,
to support [`reset_scaling!`](@ref). This function must, however, be called by the user
before assembly (to allow consistent separated assembly using different cellsets). 
```@docs
FerriteAssembly.update_scaling!
reset_scaling!
```

Additionally, [`create_threaded_scalings`](@ref) can be used to copy one scaling to each thread
when using with parallel assembly. 
```@docs
create_threaded_scalings
```

### ElementResidualScaling
```@docs
ElementResidualScaling
```
