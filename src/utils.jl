# Internal functions used for convenience

"""
    _copyto!(dest::Vector, src::Vector, inds::Vector{Int})

Internal function for faster copying of global values into the element values. 
Equivalent to `dest .= src[inds]`
"""
function _copyto!(dest::Vector, src::Vector, inds::Vector{Int})
    for (i,j) in enumerate(inds)
        dest[i] = src[j]
    end
end

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