
include("defaultvalues.jl")
include("Neumann.jl")       
include("BodyLoad.jl")
include("DofLoad.jl")

"""
    LoadHandler(dh::AbstractDofHandler)    

Create a load handler for external loads to be applied on the 
finite element simulation defined by a dofhandler `dh`. 
It can be used to apply the [`Neumann`](@ref), [`BodyLoad`](@ref),
and [`DofLoad`](@ref) contributions to the external "force"-vector, 
`fext`:

```julia
fext = zeros(ndofs(dh))
lh=LoadHandler(dh)
add!(lh, Neumann(...))  # Add Neumann boundary conditions
add!(lh, BodyLoad(...)) # Add body load
add!(lh, DofLoad(...))  # Add load to specific degrees of freedom
for t in timesteps
    fill!(fext, 0)
    apply!(fext, lh, t) # Add contributions to `fext`
    ...
end
```
"""
struct LoadHandler{DH<:Ferrite.AbstractDofHandler}
    nbcs::Dict{String,<:AbstractDomainBuffer}
    bodyloads::Dict{String,<:AbstractDomainBuffer}
    dof_loads::Vector{<:DofLoad}
    dh::DH
end
function LoadHandler(dh::Ferrite.AbstractDofHandler; threading=false)
    BT = threading ? ThreadedDomainBuffer : DomainBuffer
    return LoadHandler(Dict{String,BT}(), Dict{String,BT}(), DofLoad[], dh)
end

Ferrite.add!(lh::LoadHandler, nbc::Neumann) = (add_neumann!(lh.nbcs, nbc, lh.dh); lh)
Ferrite.add!(lh::LoadHandler, bl::BodyLoad) = (add_bodyload!(lh.bodyloads, bl, lh.dh); lh)
Ferrite.add!(lh::LoadHandler, dl::DofLoad) = (push!(lh.dof_loads, dl); lh)

# Application of boundary conditions
function Ferrite.apply!(f::Vector, lh::LoadHandler, time)
    # Abuse the Δt field, setting it to the total time.
    set_time_increment!(lh.nbcs, time)
    set_time_increment!(lh.bodyloads, time)

    worker = ReAssembler(f; fillzero=false)
    work!(worker, lh.nbcs)
    work!(worker, lh.bodyloads)
    for dof_load in lh.dof_loads
        apply_dofload!(f, dof_load, time)
    end
end
