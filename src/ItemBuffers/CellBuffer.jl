abstract type AbstractCellBuffer <: AbstractItemBuffer end

work_single!(worker, cb::AbstractCellBuffer) = work_single_cell!(worker, cb)

"""
    work_single_cell!(worker, cellbuffer)

Each worker that supports a cellbuffer should overload this function.
"""
function work_single_cell! end

mutable struct CellBuffer{T,CC,CV,DR,MT,ST,UD,UC,L} <: AbstractCellBuffer
    const dofbuffer::DofBuffer{T,DR,L}# Local dofs, dof values, residual, and stiffness
    const coords::CC                  # cellcoords (or what is required to reinit cellvalues)
    const cellvalues::CV              # Can also be named tuple of cellvalues (will be updated to MultiCellValues later)
    Δt::T                             # Time step (updated at start of assembly loop)
    cellid::Int                       # Current cell nr (updated in reinit!)
    # User defined types
    const material::MT                # User material definition (used for dispatch)
    state::ST                         # State variables re-pointed to current cell during reinit!
    old_state::ST                     # Old state variables for the cell (updated in reinit!)
    const user_data::UD               # User data for the cell (used for additional information)
    const user_cache::UC              # Cache for the cell (user type) (deepcopy for each thread)
end

"""
    CellBuffer(dofbuffer::DofBuffer, coords, cellvalues, material, state, user_data=nothing) -> CellBuffer

Create a cell buffer for an element with the local degrees of freedom described by 
`dofbuffer` (see [`DofBuffer`](@ref FerriteAssembly.DofBuffer)) and cell coordinates of the same type as `coords`.
`cellvalues` are `reinit!`ed for each cell, and the `state` is updated to the 
old cell state. `material` will be passed as-is to the element. 
The `user_data::Any` and the cache from [`allocate_cell_cache`](@ref) are available to the element via the buffer input. 

!!! note "See setup_domainbuffer"
    This constructor is normally not used, and is instead called from [`setup_domainbuffer`](@ref)

"""
function CellBuffer(dofbuffer::DofBuffer, coords, cellvalues, material, state, user_data=nothing)
    Δt = NaN 
    cellid = -1
    cache = allocate_cell_cache(material, cellvalues)
    return CellBuffer(dofbuffer, coords, cellvalues, Δt, cellid, material, state, state, user_data, cache)
end

# Required functions for a custom CellBuffer (only required internally)
# TaskLocals interface (only `create_local` required for other `AbstractCellBuffer`s) (unless gather! is req.)
function create_local(cb::CellBuffer)
    dcpy = map(deepcopy, (cb.coords, cb.cellvalues, cb.Δt, cb.cellid, cb.material, cb.state, cb.old_state))
    return CellBuffer(create_local(cb.dofbuffer), dcpy..., cb.user_data, deepcopy(cb.user_cache))
end

set_time_increment!(cb::CellBuffer, Δt) = (cb.Δt=Δt)
@inline get_dofbuffer(cb::CellBuffer) = cb.dofbuffer
@inline get_layout(cb::CellBuffer) = cb.dofbuffer.layout
@inline get_material(cb::CellBuffer) = cb.material
@inline get_values(cb::CellBuffer) = cb.cellvalues

@inline Ferrite.getcoordinates(cb::CellBuffer) = cb.coords

@inline Ferrite.cellid(cb::CellBuffer) = cb.cellid

@inline get_old_state(cb::CellBuffer) = cb.old_state

@inline get_state(cb::CellBuffer) = cb.state

@inline get_time_increment(cb::CellBuffer) = cb.Δt

@inline get_user_data(cb::CellBuffer) = cb.user_data

@inline get_user_cache(cb::CellBuffer) = cb.user_cache

"""
    FerriteAssembly.allocate_cell_cache(material, cellvalues)

This function can be overloaded for the specific material to allocate 
a cache that is stored in the `AbstractCellBuffer`. This cache can be 
used to reduce allocations. Returns `nothing` by default.
"""
allocate_cell_cache(::Any, ::Any) = nothing

"""
    reinit_buffer!(cb::CellBuffer, sim::Simulation, cellnum::Int)

Reinitialize the `cb::CellBuffer` for cell number `cellnum`.
The global degree of freedom vectors `a` (current) and `aold` are used
to update the cell degree of freedom vectors in `c`.
If the global vectors are not included in `sim`, the corresponding local vectors are set to `NaN`
The element stiffness and residual are also zeroed.
"""
function reinit_buffer!(cb::CellBuffer, sim::Simulation, cellnum::Int)
    dh = get_dofhandler(sim)
    cb.cellid = cellnum
    cb.old_state = get_old_state(sim, cellnum)
    cb.state = get_state(sim, cellnum)
    reinit_buffer!(cb.dofbuffer, sim, cellnum)
    getcoordinates!(cb.coords, dh.grid, cellnum)
    reinit!(cb.cellvalues, getcells(dh.grid, cellnum), cb.coords)
    return nothing  # Ferrite's reinit! doesn't return
end

function _replace_material_with(cb::CellBuffer, new_material)
    return setproperties(cb; material = new_material)
end

"""
    AbstractCellBufferWrapper <: AbstractCellBuffer

Supertype for cell buffers that wrap another cell buffer, e.g. to add extra information.
A subtype must define `get_wrapped(buffer)`, returning the wrapped cell buffer, to which all
accessor functions are then forwarded. A subtype can overload the accessor functions it provides
itself. Additionally, `reinit_buffer!`, `create_local`, and `_replace_material_with` must be defined.
"""
abstract type AbstractCellBufferWrapper <: AbstractCellBuffer end

function get_wrapped end

for op in (:get_dofbuffer, :get_material, :get_values, :get_time_increment, :get_state, :get_old_state,
        :get_user_data, :get_user_cache, :get_coupled_buffers, :get_layout, :(Ferrite.getcoordinates), :(Ferrite.cellid))
    @eval @inline $op(cb::AbstractCellBufferWrapper) = $op(get_wrapped(cb))
end
set_time_increment!(cb::AbstractCellBufferWrapper, Δt) = set_time_increment!(get_wrapped(cb), Δt)
