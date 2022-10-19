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

Create a cell cache for an element with `n=ndofs_per_cell(dh)`,
as well as with `materialcache`. 

    CellCache(dh::MixedDofHandler, fh::FieldHandler, materialcache=nothing)

Create a cell cache for an element with `n` degrees of freedom, 
as well as with `materialcache`. Here, `n` is the number of degrees of 
freedom for the cells in `fh.cellset`.

    CellCache(dh::MixedDofHandler, materialcache=nothing)
    CellCache(dh::MixedDofHandler, materialcache::Tuple)

Returns a tuple of `CellCache`s for each `FieldHandler` in `dh.fieldhandlers`.
If `materialcache::Tuple` is given, `materialcache[i]` is assigned to the 
`CellCache` corresponding to `dh.fieldhandlers[i]`. Otherwise, `materialcache`
is assigned to each `CellCache`
"""
function CellCache(n::Int, materialcache=nothing)
    return CellCache(zeros(n), zeros(n), zeros(n), zeros(n,n), materialcache)
end
CellCache(dh::DofHandler, args...) = CellCache(ndofs_per_cell(dh), args...)
function CellCache(dh::MixedDofHandler, fh::FieldHandler, args...)
    return CellCache(ndofs_per_cell(dh, first(fh.cellset)), args...)
end
CellCache(dh::MixedDofHandler, args...) = ntuple(i->CellCache(dh, dh.fieldhandlers[i], args...), length(dh.fieldhandlers))
function CellCache(dh::MixedDofHandler, materialcache::Tuple)
    n = length(dh.fieldhandlers)
    if n != length(materialcache)
        throw(DimensionMismatch("length(materialcache)!=length(dh.fieldhandlers)"))
    end
    return ntuple(i->CellCache(dh, dh.fieldhandlers[i], materialcache), length(dh.fieldhandlers))
end


getcellcachevars(c::CellCache) = (ae=c.anew, ae_old=c.aold, re=c.re, Ke=c.Ke, materialcache=c.materialcache)

function _copyto!(dest::Vector, src::Vector, inds::Vector{Int})
    for (i,j) in enumerate(inds)
        dest[i] = src[j]
    end
end

"""
    Ferrite.reinit!(cache::CellCache, cell::CellIterator, anew, aold)

The input is a `cache`, the `cell`, and the global degree of freedom vectors
`anew` (current) and `aold`. 
This function adds the correct element dof values to `ae` and `ae_old`, 
and resets `Ke` and `re` to zero. 
"""
function Ferrite.reinit!(cache::CellCache, cell::CellIterator, anew, aold)
    _copyto!(cache.anew, anew, celldofs(cell)) # ae_new = a_new[dofs]
    _copyto!(cache.aold, aold, celldofs(cell)) # ae_old = a_old[dofs]
    fill!(cache.Ke, 0)
    fill!(cache.re, 0)
    return cache
end