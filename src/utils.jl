# Internal functions used for convenience

"""
    _copydofs!(edofs::Vector, gdofs::Vector, inds::Vector{Int})

Internal function for faster copying of global values into the element values. 
Equivalent to `edofs .= gdofs[inds]`
    _copydofs!(edofs::Vector, gdofs::Nothing, inds::Vector{Int})

Fill `edofs` with NaN
"""
function _copydofs!(edofs::Vector, gdofs::Vector, inds::Vector{Int})
    for (i,j) in enumerate(inds)
        edofs[i] = gdofs[j]
    end
end
_copydofs!(edofs::Vector, ::Nothing, args...) = fill!(edofs, NaN)

"""
    _maketuple(t, n)

If `t` is a tuple, check that length(t)=n and return `t`. 
Otherwise, return a tuple of length `n` with `t` as every element
"""
function _maketuple(t::Tuple, n::Int)
    length(t) == n || throw(DimensionMismatch("length(t)=$(length(t))!=n=$n"))
    return t
end
_maketuple(t, n::Int) = ntuple(Returns(t), n)

"""
    _makedict(d, d_keys)

If `d::Dict` then check that it has all keys in `d_keys` (if not throw KeyError).
Otherwise, return a `Dict` with keys `d_keys` and `d` as every element
"""
_makedict(d, d_keys) = Dict(key=>d for key in d_keys)
function _makedict(d::Dict, d_keys)
    all(key->in(key, keys(d)), d_keys) || throw(KeyError("d is missing keys"))
    return d
end

"""
    intersect_nothing(a, b)

A faster intersect if `cellset::Nothing`, otherwise calls Base.intersect
"""
@inline intersect_nothing(a, ::Nothing) = a
@inline intersect_nothing(a, b) = intersect(a, b)