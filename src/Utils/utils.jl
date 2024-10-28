# Public functions for convenience 
"""
    remove_dual(x::T) where {T <: Number}
    remove_dual(x::AbstractTensor{<:Any, <:Any, T}) where {T}

Removes the dual part if `T <: ForwardDiff.Dual`, extract the value part.
Typically used when assigning state variables during differentiation calls.
"""
function remove_dual end

# Scalars
remove_dual(x::ForwardDiff.Dual) = ForwardDiff.value(x)
remove_dual(x::Number) = x

# Tensors
remove_dual(x::AbstractTensor{<:Any, <:Any, <:ForwardDiff.Dual}) = Tensors._extract_value(x)
remove_dual(x::AbstractTensor) = x

# Internal functions used for convenience

"""
    _copydofs!(edofs::Vector, gdofs::Vector, inds::Vector{Int})

Internal function for faster copying of global values into the element values. 
Equivalent to `edofs .= gdofs[inds]`

    _copydofs!(edofs::Vector, gdofs::Nothing, inds::Vector{Int})

Fill `edofs` with NaN
"""
function _copydofs!(edofs::Vector, gdofs::Vector, inds::Vector{Int})
    checkbounds(edofs, 1:length(inds))
    for (i,j) in enumerate(inds)
        gdof = gdofs[j] # checkbounds + @inbounds is slower
        @inbounds edofs[i] = gdof # checkbounds is cheap for 1:length(inds)
    end
end
_copydofs!(edofs::Vector, ::Nothing, args...) = fill!(edofs, NaN)

"""
    fast_getindex(collection)

If the output is known from the element type of the collection, this can be returned 
directly. Currently, this is implemented for `AbstractDict` with `Nothing` as type,
which takes away overhead when state variables are not used. 
"""
@inline fast_getindex(::AbstractDict{<:Any,Nothing}, key) = nothing
@inline fast_getindex(x, key) = getindex(x, key)


function intersect_cellset_sort(set::Union{AbstractSet{Int},AbstractVector{Int}}, cellset)
    intersected = set===cellset ? set : intersect(set, cellset)
    return sort!(collect(intersected))
end
function intersect_cellset_sort(set::Union{AbstractSet{FacetIndex},AbstractVector{FacetIndex}}, cellset)
    intersected_set = resize!(Vector{FacetIndex}(undef, max(length(set),length(cellset))), 0)
    for (cellnr, facetnr) in set
        cellnr ∈ cellset && push!(intersected_set, FacetIndex(cellnr, facetnr))
    end
    return sort!(intersected_set; by=first)
end

function overlaps_with_cellset(set::Union{AbstractSet,AbstractVector}, cellset)
    return any(x -> first(x) ∈ cellset, set)
end
overlaps_with_cellset(::Nothing, cellset) = !isempty(cellset)

"""
    _replace_field(container, field::Val{F}, val)

Return a modified `container`, which has the same name of the type,
but may have different parametric fields. The field with symbol `F`
is replaced by `val`, while all other fields remain the same
"""
function _replace_field(container::T, ::Val{FieldName}, val) where {T, FieldName}
    W = Base.typename(T).wrapper
    return W(map(name-> name===FieldName ? val : getproperty(container, name), fieldnames(T))...)
end