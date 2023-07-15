# This file defines iterators used for looping over a grid
# Note, most code taken from Ferrite.jl PR495
# Not exported, for convenience:
import Ferrite: UpdateFlags, AbstractGrid, AbstractDofHandler, ScalarWrapper
import Ferrite: nnodes_per_cell, cellnodes!, cellcoords!
# Just overloaded:
import Ferrite: reinit!, getnodes, getcoordinates, celldofs, cellid, celldofs!, nfaces, onboundary

@static if !isdefined(Ferrite, :FaceCache)

"""
    FaceCache(grid::Grid)
    FaceCache(dh::AbstractDofHandler)
Create a cache object with pre-allocated memory for the nodes, coordinates, and dofs of a
cell suitable for looping over faces in a grid. This struct stores the current face number,
in addition to an underlying [`CellCache`](@ref). 
The cache is updated for a new face by calling `reinit!(cache, faceindex::FaceIndex)`
**Methods with `fc::FaceCache`**
 - `reinit!(fc, faceindex)`: reinitialize the cache for face `faceindex::FaceIndex`
 - `Ferrite.faceindex(fc)`: get the `FaceIndex` of the currently cached face
 - `Ferrite.faceid(fc)`: get the current faceid (`faceindex(fc)[2]`)
 - `cellid(fc)`: get the current cellid (`faceindex(fc)[1]`)
 - `getnodes(fc)`: get the global node ids of the cell
 - `getcoordinates(fc)`: get the coordinates of the cell
 - `celldofs(fc)`: get the global dof ids of the cell
 - `reinit!(fv, fc)`: reinitialize [`FaceValues`](@ref)
 
See also [`FaceIterator`](@ref).
"""
struct FaceCache{CC<:CellCache}
    cc::CC  # const for julia>1.8 
    current_faceid::ScalarWrapper{Int} 
end
FaceCache(args...) = FaceCache(CellCache(args...), ScalarWrapper(0))

function reinit!(fc::FaceCache, face::FaceIndex)
    cellid, faceid = face
    reinit!(fc.cc, cellid)
    fc.current_faceid[] = faceid
    return nothing
end

for op = (:getnodes, :getcoordinates, :cellid, :celldofs)
    eval(quote
        function Ferrite.$op(fc::FaceCache, args...; kwargs...)
            return Ferrite.$op(fc.cc, args...; kwargs...)
        end
    end)
end
@inline faceid(fc::FaceCache) = fc.current_faceid[]
#@inline celldofs!(v::Vector, fc::FaceCache) = celldofs!(v, fc.cc)
#@inline onboundary(fc::FaceCache) = onboundary(fc.cc, faceid(fc))
#@inline faceindex(fc::FaceCache) = FaceIndex(cellid(fc), faceid(fc))
@inline function reinit!(fv::FaceValues, fc::FaceCache)
    reinit!(fv, fc.cc, faceid(fc))
end

# FaceIterator
# Leaving flags undocumented as for CellIterator
"""
    FaceIterator(gridordh::Union{Grid,AbstractDofHandler}, faceset::Set{FaceIndex})
Iterate over the faces in `faceset`. 
Create a `FaceIterator` to conveniently iterate over the faces in `faceset`. 
The elements of the iterator are [`FaceCache`](@ref)s which are properly
`reinit!`ialized. See [`FaceCache`](@ref) for more details.
Looping over a `FaceIterator`, i.e.:
```julia
for fc in FaceIterator(grid, faceset)
    # ...
end
```
is thus simply convenience for the following equivalent snippet:
```julia
fc = FaceCache(grid)
for faceindex in faceset
    reinit!(fc, faceindex)
    # ...
end
"""
struct FaceIterator{FC<:FaceCache}
    fc::FC
    set::Set{FaceIndex}
end

function FaceIterator(gridordh::Union{Grid,AbstractDofHandler}, 
                      set, flags::UpdateFlags=UpdateFlags())
    if gridordh isa MixedDofHandler
        # TODO: Keep here to maintain same settings as for CellIterator
        _check_same_celltype(gridordh.grid, set)
    end
    return FaceIterator(FaceCache(gridordh, flags), set)
end

@inline _getcache(fi::FaceIterator) = fi.fc
@inline _getset(fi::FaceIterator) = fi.set

# Iterator interface
const GridIterators{C} = FaceIterator{C}

function Base.iterate(iterator::GridIterators, state_in...)
    it = iterate(_getset(iterator), state_in...)
    it === nothing && return nothing
    cellid, state_out = it
    cache = _getcache(iterator)
    reinit!(cache, cellid)
    return (cache, state_out)
end
Base.IteratorSize(::Type{<:GridIterators}) = Base.HasLength()
Base.IteratorEltype(::Type{<:GridIterators}) = Base.HasEltype()
Base.eltype(::Type{<:GridIterators{C}}) where C = C
Base.length(iterator::GridIterators) = length(_getset(iterator))

function _check_same_celltype(grid::AbstractGrid, faceset::Set{FaceIndex})
    celltype = typeof(grid.cells[first(first(faceset))])
    if !all(typeof(grid.cells[first(face)]) == celltype for face in faceset)
        error("The cells in the faceset are not all of the same celltype.")
    end
end

end