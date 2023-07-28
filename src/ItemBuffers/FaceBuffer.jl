"""
    work_single_face!(worker, facebuffer)

Each worker that supports for a facebuffer should overload this function.
"""
function work_single_face! end

mutable struct FaceBuffer{sdim,T,FV,DR,MT,UD,UC} <: AbstractItemBuffer
    const ae_old::Vector{T}           # Old element dof values
    const ae::Vector{T}               # New element dof values
    const re::Vector{T}               # Residual/force vector 
    const Ke::Matrix{T}               # Element stiffness matrix
    const dofs::Vector{Int}           # celldofs
    const coords::Vector{Vec{sdim,T}} # cellcoords
    const facevalues::FV              # Can also be named tuple of facevalues (will be updated to MultiCellValues later)
    Δt::T                             # Time step (updated at start of assembly loop)
    cellid::Int                       # Current cell nr (updated in reinit!)
    const dofrange::DR                # dof range for each field (NamedTuple)
    # User defined types
    material::MT                      # User material definition (used for dispatch)
    const user_data::UD               # User data for the cell (used for additional information)
    const user_cache::UC              # Cache for the cell (user type) (deepcopy for each thread)
end

"""
    FaceBuffer(
        numdofs::Int, numnodes::Int, ::Val{sdim}, 
        facevalues, material, state, dofrange, user_data=nothing)

Create a facebuffer for an element with `numdofs` degrees of freedom and
`numnodes` nodes with dimension `sdim`.
`facevalues` are `reinit!`ed for each cell, and the `state` is updated to the 
old cell state. `material` will be passed as-is to the element. 
The given `dofrange::NamedTuple`, `user_data::Any`, and `cache::Any` are available to the element via the buffer input. 

!!! note "See setup_assembly"
    This constructor is normally not used, and is instead called from [`setup_assembly`](@ref)

"""
function FaceBuffer(numdofs::Int, numnodes::Int, ::Val{sdim}, facevalues, material, dofrange, user_data) where sdim
    Δt = NaN 
    cellid = -1
    cache = allocate_face_cache(material, facevalues)
    return FaceBuffer(
        zeros(numdofs), zeros(numdofs), zeros(numdofs), zeros(numdofs,numdofs), 
        zeros(Int, numdofs), zeros(Vec{sdim}, numnodes), 
        facevalues, Δt, cellid, dofrange, material, user_data, cache)
end

setup_facebuffer(ad::Bool, args...; kwargs...) = setup_facebuffer(Val(ad), args...; kwargs...)
function setup_facebuffer(::Val{false}, sdh, fv, material, dofrange, user_data)
    numdofs = ndofs_per_cell(sdh)
    numnodes = Ferrite.nnodes_per_cell(sdh)
    sdim = Val(Ferrite.getdim(sdh))
    return FaceBuffer(numdofs, numnodes, sdim, cv, material, dofrange, user_data)
end

function setup_facebuffer(::Val{true}, args...)
    error("AutoDiffBuffer not implemented for FaceBuffer")
end

# Work dispatch 
work_single!(worker, fb::FaceBuffer) = work_single_face!(worker, fb)

# Acess functions
get_ae(fb::FaceBuffer) = fb.ae
get_aeold(fb::FaceBuffer) = fb.ae_old
get_re(fb::FaceBuffer) = fb.re 
get_Ke(fb::FaceBuffer) = fb.Ke 
get_material(fb::FaceBuffer) = fb.material
get_values(fb::FaceBuffer) = fb.facevalues
get_time_increment(fb::FaceBuffer) = fb.Δt
get_user_data(fb::FaceBuffer) = fb.user_data
get_user_cache(fb::FaceBuffer) = fb.user_cache

Ferrite.celldofs(fb::FaceBuffer) = fb.dofs
Ferrite.getcoordinates(fb::FaceBuffer) = fb.coords
Ferrite.cellid(fb::FaceBuffer) = fb.cellid
Ferrite.dof_range(fb::FaceBuffer, fieldname::Symbol) = fb.dofrange[fieldname]

# Set functions 
set_time_increment!(fb::FaceBuffer, Δt) = (fb.Δt = Δt)

# TaskLocals interface
function create_local(c::FaceBuffer)
    dcpy = map(deepcopy, (c.ae_old, c.ae, c.re, c.Ke, c.dofs, c.coords, c.cellvalues, c.Δt, c.cellid, c.dofrange, c.material))
    return FaceBuffer(dcpy..., c.user_data, deepcopy(c.user_cache))
end

allocate_face_cache(::Any, ::Any) = nothing

function reinit_buffer!(fb::FaceBuffer, db::AbstractDomainBuffer, fi::FaceIndex; anew=nothing, aold=nothing)
    cellnum, facenr = fi
    dh = get_dofhandler(db)
    grid = Ferrite.get_grid(dh)
    fb.cellid = cellnum
    celldofs!(c.dofs, dh, cellnum)
    getcoordinates!(c.coords, grid, cellnum)
    reinit!(fb.facevalues, c.coords, facenr)
    _copydofs!(fb.ae,     anew, celldofs(fb)) # ae_new .= a_new[dofs]
    _copydofs!(fb.ae_old, aold, celldofs(fb)) # ae_old .= a_old[dofs]
    fill!(c.Ke, 0)
    fill!(c.re, 0)
    return nothing  # Ferrite's reinit! doesn't return 
end
