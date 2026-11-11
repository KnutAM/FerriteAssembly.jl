abstract type AbstractFacetBuffer <: AbstractItemBuffer end

work_single!(worker, fb::AbstractFacetBuffer) = work_single_facet!(worker, fb)

"""
    work_single_facet!(worker, facetbuffer)

Each worker that supports for a facetbuffer should overload this function.
"""
function work_single_facet! end


mutable struct FacetBuffer{T,CC,FV,DR,MT,UD,UC} <: AbstractFacetBuffer
    const dofbuffer::DofBuffer{T,DR}  # Local dofs, dof values, residual, and stiffness
    const coords::CC                  # cellcoords (or what is required to reinit facetvalues, in addition to facetnr)
    const facetvalues::FV              # Can also be named tuple of facetvalues (will be updated to MultiCellValues later)
    Δt::T                             # Time step (updated at start of assembly loop)
    cellid::Int                       # Current cell nr (updated in reinit!)
    # User defined types
    material::MT                      # User material definition (used for dispatch)
    const user_data::UD               # User data for the cell (used for additional information)
    const user_cache::UC              # Cache for the cell (user type) (deepcopy for each thread)
end

"""
    FacetBuffer(dofbuffer::DofBuffer, coords, facetvalues, material, user_data)

Create a facetbuffer for an element with the local degrees of freedom described by 
`dofbuffer` (see [`DofBuffer`](@ref FerriteAssembly.DofBuffer)) and cell coordinates of the same type as `coords`.
`facetvalues` are `reinit!`ed for each facet. `material` will be passed as-is to the element. 
The `user_data::Any` and the cache from [`allocate_facet_cache`](@ref) are available to the element via the buffer input. 

!!! note "See setup_domainbuffer"
    This constructor is normally not used, and is instead called from [`setup_domainbuffer`](@ref)

"""
function FacetBuffer(dofbuffer::DofBuffer, coords, facetvalues, material, user_data)
    Δt = NaN 
    cellid = -1
    cache = allocate_facet_cache(material, facetvalues)
    return FacetBuffer(dofbuffer, coords, facetvalues, Δt, cellid, material, user_data, cache)
end

# Access functions
get_dofbuffer(fb::FacetBuffer) = fb.dofbuffer
get_material(fb::FacetBuffer) = fb.material
get_values(fb::FacetBuffer) = fb.facetvalues
get_time_increment(fb::FacetBuffer) = fb.Δt
get_user_data(fb::FacetBuffer) = fb.user_data
get_user_cache(fb::FacetBuffer) = fb.user_cache

Ferrite.getcoordinates(fb::FacetBuffer) = fb.coords
Ferrite.cellid(fb::FacetBuffer) = fb.cellid

# Set functions 
set_time_increment!(fb::FacetBuffer, Δt) = (fb.Δt = Δt)

# TaskLocals interface
function create_local(c::FacetBuffer)
    dcpy = map(deepcopy, (c.coords, c.facetvalues, c.Δt, c.cellid, c.material))
    return FacetBuffer(create_local(c.dofbuffer), dcpy..., c.user_data, deepcopy(c.user_cache))
end

"""
    FerriteAssembly.allocate_facet_cache(material, cellvalues)

This function can be overloaded for the specific material to allocate 
a cache that is stored in the `FacetBuffer` and can be used to reduce 
allocations. Returns `nothing` by default.
"""
allocate_facet_cache(::Any, ::Any) = nothing

function reinit_buffer!(fb::FacetBuffer, sim::Simulation, fi::FacetIndex)
    cellnum = first(fi)
    dh = get_dofhandler(sim)
    fb.cellid = cellnum
    reinit_buffer!(fb.dofbuffer, sim, cellnum)
    reinit_values!(fb.facetvalues, fb.coords, dh.grid, fi)
    return nothing  # Ferrite's reinit! doesn't return 
end

function _replace_material_with(fb::FacetBuffer, new_material)
    return setproperties(fb; material = new_material)
end