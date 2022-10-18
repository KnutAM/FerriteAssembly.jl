struct CellCache{T,MC}
    aold::Vector{T}     # Old element dof values
    anew::Vector{T}     # New element dof values
    re::Vector{T}       # Residual/force vector 
    Ke::Matrix{T}       # Element stiffness matrix
    materialcache::MC   # Material cache
end

"""
    CellCache(n::Int, materialcache=nothing)

Create a cell cache for an element with `n` degrees of freedom,
as well as with `materialcache`

    CellCache(dh::DofHandler, materialcache=nothing)

Create a cell cache with the number of element degrees of freedom 
according to the dof handler `dh`, as well as with `materialcache`

    CellCache(dh::MixedDofHandler, fh::FieldHandler, materialcache=nothing)

Create a cell cache with the number of element degrees of freedom 
according to the mixed dof handler `dh` for its fieldhandler `fh`.
Include also the `materialcache`
"""
function CellCache(n::Int, materialcache=nothing)
    return CellCache(zeros(n), zeros(n), zeros(n), zeros(n,n), materialcache)
end
CellCache(dh::DofHandler, args...) = CellCache(ndofs_per_cell(dh), args...)
function CellCache(dh::MixedDofHandler, fh::FieldHandler, args...)
    return CellCache(ndofs_per_cell(dh, first(fh.cellset)), args...)
end

getcellcachevars(c::CellCache) = (c.anew, c.aold, c.re, c.Ke, c.materialcache)

function _copyto!(dest::Vector, src::Vector, inds::Vector{Int})
    for (i,j) in enumerate(inds)
        dest[i] = src[j]
    end
end

function Ferrite.reinit!(cache::CellCache, cell::CellIterator, anew, aold)
    _copyto!(cache.anew, anew, celldofs(cell)) # ae_new = a_new[dofs]
    _copyto!(cache.aold, aold, celldofs(cell)) # ae_old = a_old[dofs]
    fill!(cache.Ke, 0)
    fill!(cache.re, 0)
    return cache
end