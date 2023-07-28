"""
    work_single_cell!(worker, cellbuffer)

Each worker that supports a cellbuffer should overload this function.
"""
function work_single_cell! end

abstract type AbstractCellBuffer <: AbstractItemBuffer end

mutable struct CellBuffer{sdim,T,CV,DR,MT,ST,UD,UC} <: AbstractCellBuffer
    const ae_old::Vector{T}           # Old element dof values
    const ae::Vector{T}               # New element dof values
    const re::Vector{T}               # Residual/force vector 
    const Ke::Matrix{T}               # Element stiffness matrix
    const dofs::Vector{Int}           # celldofs
    const coords::Vector{Vec{sdim,T}} # cellcoords
    const cellvalues::CV              # Can also be named tuple of cellvalues (will be updated to MultiCellValues later)
    Δt::T                             # Time step (updated at start of assembly loop)
    cellid::Int                       # Current cell nr (updated in reinit!)
    const dofrange::DR                # dof range for each field (NamedTuple)
    # User defined types
    material::MT                      # User material definition (used for dispatch)
    new_state::ST                     # State variables re-pointed to current cell during reinit!
    old_state::ST                     # Old state variables for the cell (updated in reinit!)
    const user_data::UD               # User data for the cell (used for additional information)
    const user_cache::UC                   # Cache for the cell (user type) (deepcopy for each thread)
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

!!! note "See setup_assembly"
    This constructor is normally not used, and is instead called from [`setup_assembly`](@ref)

"""
function CellBuffer(numdofs::Int, numnodes::Int, ::Val{sdim}, cellvalues, material, state, dofrange, user_data=nothing) where sdim
    Δt = NaN 
    cellid = -1
    cache = allocate_cell_cache(material, cellvalues)
    return CellBuffer(
        zeros(numdofs), zeros(numdofs), zeros(numdofs), zeros(numdofs,numdofs), 
        zeros(Int, numdofs), zeros(Vec{sdim}, numnodes), 
        cellvalues, Δt, cellid, dofrange, material, state, state, user_data, cache)
end

setup_cellbuffer(ad::Bool, args...; kwargs...) = setup_cellbuffer(Val(ad), args...; kwargs...)
function setup_cellbuffer(::Val{false}, sdh, cv, material, cell_state, dofrange, user_data)
    numdofs = ndofs_per_cell(sdh)
    numnodes = Ferrite.nnodes_per_cell(sdh)
    sdim = Val(Ferrite.getdim(sdh))
    return CellBuffer(numdofs, numnodes, sdim, cv, material, cell_state, dofrange, user_data)
end

function setup_cellbuffer(::Val{true}, args...)
    return AutoDiffCellBuffer(setup_cellbuffer(Val(false), args...))
end

work_single!(worker, cb::AbstractCellBuffer) = work_single_cell!(worker, cb)

# Required functions for a custom CellBuffer (only required internally)
# TaskLocals interface (only `create_local` required for other `AbstractCellBuffer`s) (unless gather! is req.)
function create_local(cb::CellBuffer)
    dcpy = map(deepcopy, (cb.ae_old, cb.ae, cb.re, cb.Ke, cb.dofs, cb.coords, cb.cellvalues, cb.Δt, cb.cellid, cb.dofrange, cb.material, cb.new_state, cb.old_state))
    return CellBuffer(dcpy..., cb.user_data, deepcopy(cb.user_cache))
end

set_time_increment!(cb::CellBuffer, Δt) = (cb.Δt=Δt)
@inline get_Ke(cb::CellBuffer) = cb.Ke
@inline get_re(cb::CellBuffer) = cb.re
@inline get_ae(cb::CellBuffer) = cb.ae
@inline get_material(cb::CellBuffer) = cb.material
@inline get_values(cb::CellBuffer) = cb.cellvalues

"""
    Ferrite.celldofs(cb::CellBuffer)

`Ferrite.jl`'s `celldofs` function is overloaded on the `CellBuffer` to return 
the current cell's degree of freedom indices.
"""
@inline Ferrite.celldofs(cb::CellBuffer) = cb.dofs

# Convenience access functions to CellBuffer for use inside element routines
"""
    Ferrite.getcoordinates(cb::CellBuffer)

`Ferrite.jl`'s `getcoordinates` function is overloaded on the `CellBuffer` to return 
the current cell's nodal coordinates. 
"""
@inline Ferrite.getcoordinates(cb::CellBuffer) = cb.coords

"""
    FerriteAssembly.get_aeold(cb::CellBuffer)

Get the old element dof-values for the current cell 
(Filled by `NaN`s unless `aold` is passed to `doassemble!`)
"""
@inline get_aeold(cb::CellBuffer) = cb.ae_old

"""
    Ferrite.cellid(cb::CellBuffer)

Get the current cell id/nr
"""
@inline Ferrite.cellid(cb::CellBuffer) = cb.cellid

"""
    Ferrite.dof_range(cb::CellBuffer, name::Symbol)

Get the `dofrange::UnitRange{Int}` for the dofs pertaining to the field: `name`.
Same output as dof_range(dh::DofHandler, name), but fully type-stable. 
"""
@inline Ferrite.dof_range(cb::CellBuffer, name::Symbol) = cb.dofrange[name]

# Internal overload for now
Ferrite.getfieldnames(cb::CellBuffer) = keys(cb.dofrange)

"""
    get_old_state(cb::CellBuffer)

Get the state variables for the cell from the previous time step. 

!!! note
    If no `old_state` keyword is passed to `doassemble!`, this variable 
    will not be updated for the given cell, and typically contains the 
    initial cell state. 
"""
@inline get_old_state(cb::CellBuffer) = cb.old_state

@inline get_new_state(cb::CellBuffer) = cb.new_state

"""
    get_time_increment(cb::CellBuffer)

Get the time increment, `Δt`, that was passed to `doassemble` (defaults to `NaN`)
"""
@inline get_time_increment(cb::CellBuffer) = cb.Δt

"""
    FerriteAssembly.get_user_data(cb::CellBuffer)

Get the user specified `user_data` given to `setup_assembly` or `AssemblyDomain`
"""
@inline get_user_data(cb::CellBuffer) = cb.user_data

"""
    FerriteAssembly.get_user_cache(cb::CellBuffer)

Get the user-specified `cache` created by [`allocate_cell_cache`](@ref)
"""
@inline get_user_cache(cb::CellBuffer) = cb.user_cache

get_cache(::CellBuffer) = error("get_cache discontinued, user get_user_cache instead")

"""
    FerriteAssembly.allocate_cell_cache(material, cellvalues)

This function can be overloaded for the specific material to allocate 
a cache that is stored in the `CellBuffer` which can be used to reduce 
allocations in the element routine. Returns `nothing` by default.
"""
allocate_cell_cache(::Any, ::Any) = nothing

"""
    reinit_buffer!(cb::CellBuffer, dh::AbstractDofHandler, cellnum::Int; anew=nothing, aold=nothing)

Reinitialize the `cb::CellBuffer` for cell number `cellnum`.
The global degree of freedom vectors `anew` (current) and `aold` are used
to update the cell degree of freedom vectors in `c`.
If the global vectors are instead `::Nothing`, the corresponding cell values are set to `NaN`
The element stiffness, `cb.Ke`, and residual, `cb.re`, are also zeroed. 
"""
function reinit_buffer!(cb::CellBuffer, db::AbstractDomainBuffer, cellnum::Int; anew=nothing, aold=nothing)
    dh = get_dofhandler(db)
    grid = dh.grid
    cb.cellid = cellnum
    cb.old_state = get_old_state(db, cellnum)
    cb.new_state = get_new_state(db, cellnum)
    celldofs!(cb.dofs, dh, cellnum)
    getcoordinates!(cb.coords, grid, cellnum)
    reinit!(cb.cellvalues, cb.coords)
    _copydofs!(cb.ae,     anew, cb.dofs) # ae_new .= a_new[dofs]
    _copydofs!(cb.ae_old, aold, cb.dofs) # ae_old .= a_old[dofs]
    fill!(cb.Ke, 0)
    fill!(cb.re, 0)
    return nothing  # Ferrite's reinit! doesn't return 
end
