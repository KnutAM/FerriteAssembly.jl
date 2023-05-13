```@meta
CurrentModule = FerriteAssembly
```

## CellBuffer
Variables that are used and modified for each cell of a certain type, 
but that don't belong to a specific cell, are collected in a `CellBuffer`.

Note that construction of `CellBuffer` happens automatically when calling
[`setup_assembly`](@ref), and to include `cache` or `user_data`, these 
should be passed to `setup_assembly` (or [`AssemblyDomain`](@ref)).

### Access functions
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