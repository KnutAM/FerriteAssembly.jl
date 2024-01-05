```@meta
CurrentModule = FerriteAssembly
```

# AbstractItemBuffer
Depending on the domain that is worked over, different item buffers are available, 
e.g. `CellBuffer` and `FaceBuffer`. These are set up during call to `setup_domainbuffer`.

For each item (cell or face), the values in the buffer are updated to the current item,
and can be accessed with the following functions:

```@docs
get_aeold
get_old_state
get_time_increment
Ferrite.dof_range(::AbstractItemBuffer, ::Symbol)
Ferrite.getcoordinates(::AbstractItemBuffer)
Ferrite.celldofs(::AbstractItemBuffer)
Ferrite.cellid(::AbstractItemBuffer)
get_user_data
get_user_cache
```

## AbstractCellBuffer
Two types of cell buffers are provided, the regular `CellBuffer` and a wrapper `AutoDiffCellBuffer` that speeds up the automatic differentiation. By using the above access functions, these behave identical.

To allocate the user cache, overload `allocate_cell_cache`:
```@docs
allocate_cell_cache(::Any, ::Any)
```

## FaceBuffer
The `FaceBuffer` is similar to the `CellBuffer`, except that it has `FaceValues` instead of `CellValues`
and that state variables are not supported. To allocate user cache, overload `allocate_face_cache`:
```@docs
allocate_face_cache
```