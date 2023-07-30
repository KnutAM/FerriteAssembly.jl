
include("FerritePR495.jl")  # FaceIterator
include("defaultvalues.jl")
include("Neumann.jl")       
include("BodyLoad.jl")

"""
    LoadHandler(dh::AbstractDofHandler)    

The load handler for external loads on a dofhandler `dh`. 
It can be used to apply the Neumann and body load  
contributions to the external "force"-vector, `fext`:

```julia
fext = zeros(ndofs(dh))
lh=LoadHandler(dh)
add!(lh, Neumann(...))  # Add Neumann boundary conditions
add!(lh, BodyLoad(...)) # Add body load
for t in timesteps
    fill!(fext, 0)
    apply!(fext, lh, t) # Add contributions to `fext`
    ...
end
```
"""
struct LoadHandler{DH<:AbstractDofHandler}
    nbcs::Dict{String,<:AbstractDomainBuffer}
    bodyloads::Dict{String,<:AbstractDomainBuffer}
    dh::DH
end
function LoadHandler(dh::AbstractDofHandler; threading=false)
    BT = threading ? ThreadedDomainBuffer : DomainBuffer
    return LoadHandler(Dict{String,BT}(), Dict{String,BT}(), dh)
end

Ferrite.add!(lh::LoadHandler, nbc::Neumann) = add_neumann!(lh.nbcs, nbc, lh.dh)
Ferrite.add!(lh::LoadHandler, bl::BodyLoad) = add_bodyload!(lh.bodyloads, bl, lh.dh)

# Application of boundary conditions
function Ferrite.apply!(f::Vector, lh::LoadHandler, time)
    # Abuse the Δt field, setting it to the total time.
    set_time_increment!(lh.nbcs, time)
    set_time_increment!(lh.bodyloads, time)

    worker = ReAssembler(f; fillzero=false)
    work!(worker, lh.nbcs)
    work!(worker, lh.bodyloads)
end
