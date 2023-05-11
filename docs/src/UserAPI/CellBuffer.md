```@meta
CurrentModule = FerriteAssembly
```

## CellBuffer
Variables that are used and modified for each cell of a certain type, 
but that don't belong to a specific cell, are collected in a `CellBuffer`.

### Access functions
The following access functions can be used to extract information from the 
`CellBuffer`
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