```@meta
CurrentModule = FerriteAssembly
```

# API
Methods that should be overloaded and exported functions are described on this page.

## Setting up and doing the assembly
```@docs
setup_assembly
AssemblyDomain
doassemble!
```

## Element routines
One of the element methods should be overloaded for a given combination of `cellvalues`
and `material`. 
Note that the `cellvalues` are already `reinit!`:ed when passed to the element routines. 
```@docs
FerriteAssembly.element_routine!
FerriteAssembly.element_residual!
```

## State variables
The initial state variables may vary depending on the position in the grid.
Furthermore, the datastructure depends on the type of dof handler, so
a convenience function exists that creates the correct variable structure. 
```@docs
FerriteAssembly.create_cell_state
update_states!
FerriteAssembly.create_states
```

## CellBuffer
Variables that are used and modified for each cell of a certain type, 
but that don't belong to a specific cell, are collected in a `CellBuffer`.

### Access functions
The following access functions can be used to extract information from the 
`CellBuffer`
```@docs
get_state_old
get_aeold
get_time_increment
Ferrite.dof_range
Ferrite.getcoordinates
Ferrite.celldofs
Ferrite.cellid
get_user_data
get_cache
```

## Residual scaling
There are many options for how to scale the residual in finite element simulations.
This package does not intend to implement many different options, but does give the 
user the option to calculate scaling contributions from each cell, which may be useful.
By defining a `scaling` that is passed to [`setup_assembly`](@ref), it can be updated 
based on the output from each cell. 

One type of scaling, [`ElementResidualScaling`](@ref), is included as described below.
Its code can be used as a template for how to include custom scaling that works on 
the element level.

In general, a scaling must support the following functions
```@docs
FerriteAssembly.update_scaling!
reset_scaling!
add_to_scaling!
```

### ElementResidualScaling
```@docs
ElementResidualScaling
```
