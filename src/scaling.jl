"""
    update_scaling!(scaling, re, cellbuffer)

This function should add the contribution from the element residual vector `re` to the
scaling factors in `scaling`.
"""
function update_scaling! end

"""
    reset_scaling!(scaling)

This function should reset the scaling factors, such that the values don't accumulate if used 
in multiple iterations/time steps. 
"""
function reset_scaling! end

"""
    add_to_scaling!(base_scaling, scaling)

Add the scaling values from `scaling` to `base_scaling`. This function is used after 
threaded assembly to put all the scaling values together into the object initially sent
to `setup_assembly`
"""
function add_to_scaling! end


# Default: nothing = no scaling 
update_scaling!(::Nothing, args...) = nothing
reset_scaling!(::Nothing, args...) = nothing 
add_to_scaling!(::Nothing, args...) = nothing

"""
    ElementResidualScaling(dh::AbstractDofHandler, p=Val(2), T=Float64)

Create tolerance scaling based on the sum of the p-norm of each cell's residual vector, 
separately for each field. I.e. pseudo-code for field `:u` the scaling `factor::T` is
```
for each cell
    element_routine!(Ke, re, args...) # Calculate re
    factor += sum(abs.(re[dof_range(dh, :u)]).^p)^(1/p)
```
Note that `p=Val(2)` is a special case that is supported, for different values 
it should be given as a `Real`. `p=2` is equivalent to `Val(2)`, but is less efficient. 
"""
struct ElementResidualScaling{T,P}
    factors::Dict{Symbol,T}
    p::P
end
function ElementResidualScaling(dh, p=Val(2), T=Float64)
    return ElementResidualScaling(Dict(key=>zero(T) for key in Ferrite.getfieldnames(dh)), p)
end

function update_scaling!(s::ElementResidualScaling, re, buffer, args...)
    p = s.p
    pinv = 1/p 
    for fieldname in keys(s.factors) 
        if fieldname in Ferrite.getfieldnames(buffer)
            s.factors[fieldname] += (sum(i->abs(re[i])^p, dof_range(buffer, fieldname)))^pinv
        end
    end
end
function update_scaling!(s::ElementResidualScaling{<:Any,Val{2}}, re, buffer, args...)
    for fieldname in keys(s.factors) 
        if fieldname in Ferrite.getfieldnames(buffer)
            s.factors[fieldname] += sqrt(sum(i->re[i]^2, dof_range(buffer, fieldname)))
        end
    end
end

function reset_scaling!(s::ElementResidualScaling{T}, args...) where T
    for fieldname in keys(s.factors)
        s.factors[fieldname] = zero(T)
    end
end

function add_to_scaling!(base_scaling::T, scaling::T) where {T<:ElementResidualScaling}
    for fieldname in keys(base_scaling.factors)
        base_scaling.factors[fieldname] += scaling.factors[fieldname]
    end
end