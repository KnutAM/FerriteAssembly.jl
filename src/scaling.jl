create_threaded_scalings(scaling; nthreads=Threads.nthreads()) = [deepcopy(scaling) for _ in 1:nthreads]

# Default, no scaling at all
update_scaling!(::Nothing, args...) = nothing
reset_scaling!(::Nothing, args...) = nothing 

"""
    ElementResidualScaling(dh::AbstractDofHandler, p, T=Float64)

Create tolerance scaling based on the sum of each cell's norm(re, p), for each field separately.
I.e. 
```
factor[fieldname] = 0.0
for each cell
    calculate re 
    factor[fieldname] += outerfun(sum(innerfun.(re[dof_range(dh, fieldname)])))
```
"""
struct ElementResidualScaling{T,P}
    factors::Dict{Symbol,T}
    p::P
end
function ElementResidualScaling(dh, p=2, T=Float64)
    return ElementResidualScaling(Dict(key=>zero(T) for key in Ferrite.getfieldnames(dh)), p)
end

function update_scaling!(s::ElementResidualScaling, re, dh_fh, args...)
    p = s.p
    pinv = 1/p 
    for fieldname in keys(s.factors)
        s.factors[fieldname] += (sum(i->abs(re[i])^p, dof_range(dh_fh, fieldname)))^pinv
    end
end

function reset_scaling!(s::ElementResidualScaling{T}, args...) where T
    for fieldname in keys(s.factors)
        s.factors[fieldname] = zero(T)
    end
end