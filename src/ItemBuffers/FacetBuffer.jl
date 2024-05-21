abstract type AbstractFacetBuffer <: AbstractItemBuffer end

work_single!(worker, fb::AbstractFacetBuffer) = work_single_facet!(worker, fb)

"""
    work_single_facet!(worker, facetbuffer)

Each worker that supports for a facetbuffer should overload this function.
"""
function work_single_facet! end


mutable struct FacetBuffer{T,CC,FV,DR,MT,UD,UC} <: AbstractFacetBuffer
    const ae_old::Vector{T}           # Old element dof values
    const ae::Vector{T}               # New element dof values
    const re::Vector{T}               # Residual/force vector 
    const Ke::Matrix{T}               # Element stiffness matrix
    const dofs::Vector{Int}           # celldofs
    const coords::CC                  # cellcoords (or what is required to reinit facetvalues, in addition to facetnr)
    const facetvalues::FV              # Can also be named tuple of facetvalues (will be updated to MultiCellValues later)
    Δt::T                             # Time step (updated at start of assembly loop)
    cellid::Int                       # Current cell nr (updated in reinit!)
    const dofrange::DR                # dof range for each field (NamedTuple)
    # User defined types
    material::MT                      # User material definition (used for dispatch)
    const user_data::UD               # User data for the cell (used for additional information)
    const user_cache::UC              # Cache for the cell (user type) (deepcopy for each thread)
end

"""
    FacetBuffer(
        numdofs::Int, numnodes::Int, ::Val{sdim}, 
        facetvalues, material, state, dofrange, user_data=nothing)

Create a facetbuffer for an element with `numdofs` degrees of freedom and
`numnodes` nodes with dimension `sdim`.
`facetvalues` are `reinit!`ed for each cell, and the `state` is updated to the 
old cell state. `material` will be passed as-is to the element. 
The given `dofrange::NamedTuple`, `user_data::Any`, and `cache::Any` are available to the element via the buffer input. 

!!! note "See setup_domainbuffer"
    This constructor is normally not used, and is instead called from [`setup_domainbuffer`](@ref)

"""
function FacetBuffer(numdofs::Int, coords, facetvalues, material, dofrange, user_data)
    Δt = NaN 
    cellid = -1
    cache = allocate_facet_cache(material, facetvalues)
    return FacetBuffer(
        zeros(numdofs), zeros(numdofs), zeros(numdofs), zeros(numdofs,numdofs), 
        zeros(Int, numdofs), coords, 
        facetvalues, Δt, cellid, dofrange, material, user_data, cache)
end

setup_facetbuffer(ad::Bool, args...; kwargs...) = setup_facetbuffer(Val(ad), args...; kwargs...)
function setup_facetbuffer(::Val{false}, sdh, fv, material, dofrange, user_data)
    numdofs = ndofs_per_cell(sdh)
    coords = getcoordinates(_getgrid(sdh), first(_getcellset(sdh)))
    return FacetBuffer(numdofs, coords, fv, material, dofrange, user_data)
end

function setup_facetbuffer(::Val{true}, args...)
    error("AutoDiffBuffer not implemented for FacetBuffer")
end

# Access functions
get_ae(fb::FacetBuffer) = fb.ae
get_aeold(fb::FacetBuffer) = fb.ae_old
get_re(fb::FacetBuffer) = fb.re 
get_Ke(fb::FacetBuffer) = fb.Ke 
get_material(fb::FacetBuffer) = fb.material
get_values(fb::FacetBuffer) = fb.facetvalues
get_time_increment(fb::FacetBuffer) = fb.Δt
get_user_data(fb::FacetBuffer) = fb.user_data
get_user_cache(fb::FacetBuffer) = fb.user_cache

Ferrite.celldofs(fb::FacetBuffer) = fb.dofs
Ferrite.getcoordinates(fb::FacetBuffer) = fb.coords
Ferrite.cellid(fb::FacetBuffer) = fb.cellid
Ferrite.dof_range(fb::FacetBuffer, fieldname::Symbol) = fb.dofrange[fieldname]

# Set functions 
set_time_increment!(fb::FacetBuffer, Δt) = (fb.Δt = Δt)

# Internal overload 
Ferrite.getfieldnames(fb::FacetBuffer) = keys(fb.dofrange)

# TaskLocals interface
function create_local(c::FacetBuffer)
    dcpy = map(deepcopy, (c.ae_old, c.ae, c.re, c.Ke, c.dofs, c.coords, c.facetvalues, c.Δt, c.cellid, c.dofrange, c.material))
    return FacetBuffer(dcpy..., c.user_data, deepcopy(c.user_cache))
end

"""
    FerriteAssembly.allocate_facet_cache(material, cellvalues)

This function can be overloaded for the specific material to allocate 
a cache that is stored in the `FacetBuffer` and can be used to reduce 
allocations. Returns `nothing` by default.
"""
allocate_facet_cache(::Any, ::Any) = nothing

function reinit_buffer!(fb::FacetBuffer, db::AbstractDomainBuffer, fi::FacetIndex; a=nothing, aold=nothing)
    cellnum, facetnr = fi
    dh = get_dofhandler(db)
    fb.cellid = cellnum
    celldofs!(fb.dofs, dh, cellnum)
    getcoordinates!(fb.coords, dh.grid, cellnum)
    reinit!(fb.facetvalues, fb.coords, facetnr)
    _copydofs!(fb.ae,     a,    celldofs(fb)) # ae_new .= a_new[dofs]
    _copydofs!(fb.ae_old, aold, celldofs(fb)) # ae_old .= a_old[dofs]
    fill!(fb.Ke, 0)
    fill!(fb.re, 0)
    return nothing  # Ferrite's reinit! doesn't return 
end

function _replace_material_with(fb::FacetBuffer, new_material)
    return _replace_field(fb, Val(:material), new_material)
end