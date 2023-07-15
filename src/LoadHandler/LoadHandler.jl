
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
    nbcs::Vector{NeumannData}
    bodyloads::Vector{BodyLoadData}
    dh::DH
end
LoadHandler(dh::AbstractDofHandler) = LoadHandler(NeumannData[], BodyLoadData[], dh)

Ferrite.add!(lh::LoadHandler, nbc::Neumann) = add_neumann!(lh.nbcs, nbc, lh.dh)
Ferrite.add!(lh::LoadHandler, bl::BodyLoad) = add_bodyload!(lh.bodyloads, bl, lh.dh)

# Application of boundary conditions
function Ferrite.apply!(f::Vector, lh::LoadHandler, time)
    foreach(nbc->apply_neumann!(f,nbc,lh.dh,time), lh.nbcs)
    foreach(bl->apply_bodyload!(f,bl,lh.dh,time), lh.bodyloads)
end
