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
AutoDiffCellBuffer
getCellBuffer
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