abstract type AbstractCellBuffer <: AbstractItemBuffer end

work_single!(worker, cb::AbstractCellBuffer) = work_single_cell!(worker, cb)

"""
    work_single_cell!(worker, cellbuffer)

Each worker that supports a cellbuffer should overload this function.
"""
function work_single_cell! end

mutable struct CellBuffer{T,CC,CV,DR,MT,ST,UD,UC,CB} <: AbstractCellBuffer
    const ae_old::Vector{T}           # Old element dof values
    const ae::Vector{T}               # Current element dof values
    const re::Vector{T}               # Residual/force vector 
    const Ke::Matrix{T}               # Element stiffness matrix
    const dofs::Vector{Int}           # celldofs
    const coords::CC                  # cellcoords (or what is required to reinit cellvalues)
    const cellvalues::CV              # Can also be named tuple of cellvalues (will be updated to MultiCellValues later)
    Δt::T                             # Time step (updated at start of assembly loop)
    cellid::Int                       # Current cell nr (updated in reinit!)
    const dofrange::DR                # dof range for each field (NamedTuple)
    # User defined types
    const material::MT                # User material definition (used for dispatch)
    state::ST                         # State variables re-pointed to current cell during reinit!
    old_state::ST                     # Old state variables for the cell (updated in reinit!)
    const user_data::UD               # User data for the cell (used for additional information)
    const user_cache::UC              # Cache for the cell (user type) (deepcopy for each thread)
    const coupled_buffers::CB         # nothing or NamedTuple with staggered coupled `CellBuffer`s. 
end

"""
    CellBuffer(
        numdofs::Int, numnodes::Int, ::Val{sdim}, 
        cellvalues, material, state, dofrange, user_data=nothing) -> CellBuffer

Create a cell cache for an element with `numdofs` degrees of freedom and
`numnodes` nodes with dimension `sdim`.
`cellvalues` are `reinit!`ed for each cell, and the `state` is updated to the 
old cell state. `material` will be passed as-is to the element. 
The given `dofrange::NamedTuple`, `user_data::Any`, and `cache::Any` are available to the element via the buffer input. 

!!! note "See setup_domainbuffer"
    This constructor is normally not used, and is instead called from [`setup_domainbuffer`](@ref)

"""
function CellBuffer(numdofs::Int, coords, cellvalues, material, state, dofrange, user_data=nothing)
    Δt = NaN 
    cellid = -1
    cache = allocate_cell_cache(material, cellvalues)
    return CellBuffer(
        zeros(numdofs), zeros(numdofs), zeros(numdofs), zeros(numdofs,numdofs), 
        zeros(Int, numdofs), coords, 
        cellvalues, Δt, cellid, dofrange, material, state, state, user_data, cache, nothing)
end

setup_cellbuffer(ad::Bool, args...; kwargs...) = setup_cellbuffer(Val(ad), args...; kwargs...)
function setup_cellbuffer(::Val{false}, sdh, cv, material, cell_state, dofrange, user_data)
    numdofs = ndofs_per_cell(sdh)
    coords = getcoordinates(_getgrid(sdh), first(_getcellset(sdh)))
    return CellBuffer(numdofs, coords, cv, material, cell_state, dofrange, user_data)
end

function couple_cellbuffers(;kwargs...)
    return map((k, cb) -> couple_cellbuffers(cb; kwargs...), pairs(kwargs))
end
function couple_cellbuffers(cb::CB; kwargs...) where {CB <: CellBuffer}
    # TODO: Improve type-stability (perhaps not critical?)
    # Possible to pass the key as a Val type if required...
    coupled_buffers = NamedTuple(k => v for (k, v) in kwargs if v !== cb)
    return  _replace_field(cb, Val(:coupled_buffers), coupled_buffers)
end

struct_to_namedtuple(x::T) where T = NamedTuple{fieldnames(T)}(tuple((getproperty(x,k) for k in fieldnames(T))...));
function setup_cellbuffer(::Val{true}, args...)
    return AutoDiffCellBuffer(setup_cellbuffer(Val(false), args...))
end

# Required functions for a custom CellBuffer (only required internally)
# TaskLocals interface (only `create_local` required for other `AbstractCellBuffer`s) (unless gather! is req.)
function create_local(cb::CellBuffer)
    dcpy = map(deepcopy, (cb.ae_old, cb.ae, cb.re, cb.Ke, cb.dofs, cb.coords, cb.cellvalues, cb.Δt, cb.cellid, cb.dofrange, cb.material, cb.state, cb.old_state))
    return CellBuffer(dcpy..., cb.user_data, deepcopy(cb.user_cache), create_local(cb.coupled_buffers))
end

set_time_increment!(cb::CellBuffer, Δt) = (cb.Δt=Δt)
@inline get_Ke(cb::CellBuffer) = cb.Ke
@inline get_re(cb::CellBuffer) = cb.re
@inline get_ae(cb::CellBuffer) = cb.ae
@inline get_material(cb::CellBuffer) = cb.material
@inline get_values(cb::CellBuffer) = cb.cellvalues

@inline Ferrite.celldofs(cb::CellBuffer) = cb.dofs

@inline Ferrite.getcoordinates(cb::CellBuffer) = cb.coords

@inline get_aeold(cb::CellBuffer) = cb.ae_old

@inline Ferrite.cellid(cb::CellBuffer) = cb.cellid

@inline Ferrite.dof_range(cb::CellBuffer, name::Symbol) = cb.dofrange[name]

# Internal overload for now
Ferrite.getfieldnames(cb::CellBuffer) = keys(cb.dofrange)

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
    reinit_buffer!(cb::CellBuffer, db::AbstractDomainBuffer, cellnum::Int; a=nothing, aold=nothing)

Reinitialize the `cb::CellBuffer` for cell number `cellnum`.
The global degree of freedom vectors `a` (current) and `aold` are used
to update the cell degree of freedom vectors in `c`.
If the global vectors are instead `::Nothing`, the corresponding cell values are set to `NaN`
The element stiffness, `cb.Ke`, and residual, `cb.re`, are also zeroed. 
"""
function reinit_buffer!(cb::CellBuffer, db::AbstractDomainBuffer, cellnum::Int; a=nothing, aold=nothing)
    dh = get_dofhandler(db)
    grid = dh.grid
    cb.cellid = cellnum
    cb.old_state = get_old_state(db, cellnum)
    cb.state = get_state(db, cellnum)
    celldofs!(cb.dofs, dh, cellnum)
    getcoordinates!(cb.coords, grid, cellnum)
    reinit!(cb.cellvalues, getcells(grid, cellnum), cb.coords)
    _copydofs!(cb.ae,     a, cb.dofs) # ae_new .= a_new[dofs]
    _copydofs!(cb.ae_old, aold, cb.dofs) # ae_old .= a_old[dofs]
    fill!(cb.Ke, 0)
    fill!(cb.re, 0)
    return nothing  # Ferrite's reinit! doesn't return 
end

function _replace_material_with(cb::CellBuffer, new_material)
    return setproperties(cb; material = new_material)
end
