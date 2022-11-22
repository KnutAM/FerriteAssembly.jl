create_threaded_scalings(scaling; nthreads=Threads.nthreads()) = [deepcopy(scaling) for _ in 1:nthreads]

# Default, no scaling at all
update_scaling!(::Nothing, args...) = nothing
reset_scaling!(::Nothing, args...) = nothing 

"""
    ElementResidualScaling(dh::AbstractDofHandler, T=Float64)

Create tolerance scaling based on the sum of each cell's norm(re), for each field separately.
I.e. 
```
factor[fieldname] = 0.0
for each cell
    calculate re 
    factor[fieldname] += norm(re[dof_range(dh, fieldname)])
```
"""
struct ElementResidualScaling{T}
    factors::Dict{Symbol,T}
end
ElementResidualScaling(dh, T=Float64) = ElementResidualScaling(Dict(key=>zero(T) for key in Ferrite.getfieldnames(dh)))

function update_scaling!(s::ElementResidualScaling, re, dh_fh, args...)
    for fieldname in keys(s.factors)
        s.factors[fieldname] += sqrt(sum(i->re[i]^2, dof_range(dh_fh, fieldname)))
    end
end

function reset_scaling!(s::ElementResidualScaling{T}, args...) where T
    for fieldname in keys(s.factors)
        s.factors[fieldname] = zero(T)
    end
end