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
    intersect_nothing(a, b)

A faster intersect if `cellset::Nothing`, otherwise calls Base.intersect
"""
@inline intersect_nothing(a, ::Nothing) = a
@inline intersect_nothing(a, b) = intersect(a, b)