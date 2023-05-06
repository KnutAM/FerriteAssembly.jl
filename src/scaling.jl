"""
    update_scaling!(scaling, re, dh_fh, cellbuffer)

This function should add the contribution from the element residual vector `re` to the
scaling factors in `scaling`. 
The cellbuffer is included, which contains extra information if needed. 
"""
function update_scaling! end

"""
    reset_scaling!(scaling)

This function should reset the scaling factors, such that the values don't accumulate if used 
in multiple iterations/time steps. 
"""
function reset_scaling! end

"""
    create_threaded_scalings(scaling; nthreads=Threads.nthreads()) -> Vector

Makes `nthreads` deepcopies of `scaling`, such that it can be used for threaded assembly.
Note: If creating your own scaling, it might make sense to define 
`Base.sum(::Vector{<:MyScalingType})` to sum scalings after threaded assembly. 
"""
create_threaded_scalings(scaling; nthreads=Threads.nthreads()) = [deepcopy(scaling) for _ in 1:nthreads]

# Default: nothing = no scaling 
update_scaling!(::Nothing, args...) = nothing
reset_scaling!(::Nothing, args...) = nothing 

"""
    ElementResidualScaling(dh::AbstractDofHandler, p=2, T=Float64)

Create tolerance scaling based on the sum of the p-norm of each cell's residual vector, 
separately for each field. I.e. pseudo-code for field `:u` the scaling `factor::T` is
```
for each cell
    element_routine!(Ke, re, args...) # Calculate re
    factor += sum(abs.(re[dof_range(dh, :u)]).^p)^(1/p)
```
"""
struct ElementResidualScaling{T,P}
    factors::Dict{Symbol,T}
    p::P
end
function ElementResidualScaling(dh, p=2, T=Float64)
    return ElementResidualScaling(Dict(key=>zero(T) for key in Ferrite.getfieldnames(dh)), p)
end

function update_scaling!(s::ElementResidualScaling, re, buffer, args...)
    p = s.p
    pinv = 1/p 
    for fieldname in keys(s.factors)
        if haskey(dof_range(buffer), fieldname)
            s.factors[fieldname] += (sum(i->abs(re[i])^p, dof_range(buffer, fieldname)))^pinv
        end
    end
end

function reset_scaling!(s::ElementResidualScaling{T}, args...) where T
    for fieldname in keys(s.factors)
        s.factors[fieldname] = zero(T)
    end
end

Base.sum(scaling::ElementResidualScaling) = scaling     # ok as Base.sum(v::Number)=v is defined...
function Base.sum(scalings::Vector{<:ElementResidualScaling})
    _getfactor(s::ElementResidualScaling, fieldname) = s.factors[fieldname]

    factors = Dict(key=>sum(s->_getfactor(s,key), scalings) for key in keys(first(scalings).factors))
    return ElementResidualScaling(factors, first(scalings).p)
end
