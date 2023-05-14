```@meta
CurrentModule = FerriteAssembly
```

# CellBuffer
A `CellBuffer` contains variables that are used and modified for each cell of a certain type.
For each cell, it is `reinit!`:ed such that values such as cell coordinates, old state variables, 
etc. are updated to the current cell. 

## Construction
Construction of `CellBuffer` happens automatically when calling
[`setup_assembly`](@ref). 
To extra data can be included via the `user_data` keyword when calling
`setup_assembly` or [`AssemblyDomain`](@ref).

To allocate cache for the element routine, overload `allocate_cell_cache`:
```@docs
allocate_cell_cache(::Any, ::Any)
```

## Access functions
The following access functions can be used to extract information from 
an `AbstractCellBuffer`.
```@docs
get_old_state
get_aeold
get_time_increment
Ferrite.dof_range
Ferrite.getcoordinates
Ferrite.celldofs
Ferrite.cellid
get_user_data
get_cache
```