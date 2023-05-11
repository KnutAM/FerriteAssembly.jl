abstract type AbstractCellBuffer end

mutable struct CellBuffer{sdim,T,CV,DR,MT,ST,UD,CC} <: AbstractCellBuffer
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
    const material::MT                # User material definition (used for dispatch)
    old_state::ST                     # Old state variables for the cell (updated in reinit!)
    const user_data::UD               # User data for the cell (used for additional information)
    const cache::CC                   # Cache for the cell (user type) (deepcopy for each thread)
end

"""
    CellBuffer(
        numdofs::Int, numnodes::Int, ::Val{sdim}, 
        cellvalues, material, cell_load=nothing, cache=nothing) -> CellBuffer

Create a cell cache for an element with `numdofs` degrees of freedom and
`numnodes` nodes with dimension `sdim`.
`cellvalues` are `reinit!`ed for each cell, and the `state` is updated to the 
old cell state. `material` will be passed as-is to the element. 
The given `dofrange::NamedTuple`, `user_data::Any`, and `cache::Any` are available to the element via the buffer input. 

!!! note "See setup_assembly"
    This constructor is normally not used, and is instead called from [`setup_assembly`](@ref)

"""
function CellBuffer(numdofs::Int, numnodes::Int, ::Val{sdim}, cellvalues, material, state, dofrange, user_data=nothing, cache=nothing) where sdim
    Δt = NaN 
    cellid = -1
    return CellBuffer(
        zeros(numdofs), zeros(numdofs), zeros(numdofs), zeros(numdofs,numdofs), 
        zeros(Int, numdofs), zeros(Vec{sdim}, numnodes), 
        cellvalues, Δt, cellid, dofrange, material, state, user_data, cache)
end

setup_cellbuffer(ad::Bool, args...; kwargs...) = setup_cellbuffer(Val(ad), args...; kwargs...)
function setup_cellbuffer(::Val{false}, sdh, cv, material, cell_state, dofrange, user_data, cache)
    numdofs = ndofs_per_cell(sdh)
    numnodes = Ferrite.nnodes_per_cell(sdh)
    sdim = Val(Ferrite.getdim(sdh))
    return CellBuffer(numdofs, numnodes, sdim, cv, material, cell_state, dofrange, user_data, cache)
end

function setup_cellbuffer(::Val{true}, args...)
    return AutoDiffCellBuffer(setup_cellbuffer(Val(false), args...))
end

# Required functions for a custom CellBuffer (only required internally)
# TaskLocals interface (only `create_local` required for other `AbstractCellBuffer`s) (unless gather! is req.)
function create_local(c::CellBuffer)
    dcpy = (deepcopy(f) for f in (c.ae_old, c.ae, c.re, c.Ke, c.dofs, c.coords, c.cellvalues, c.Δt, c.cellid, c.dofrange, c.material, c.old_state))
    return CellBuffer(dcpy..., c.user_data, deepcopy(c.cache))
end
scatter!(task::AbstractCellBuffer, base::AbstractCellBuffer) = update_time!(task, get_time_increment(base))
gather!(::AbstractCellBuffer, ::AbstractCellBuffer) = nothing

update_time!(c::CellBuffer, Δt) = setfield!(c, :Δt, Δt)
@inline get_Ke(c::CellBuffer) = c.Ke
@inline get_re(c::CellBuffer) = c.re
@inline get_ae(c::CellBuffer) = c.ae
@inline get_material(c::CellBuffer) = c.material
@inline get_cellvalues(c::CellBuffer) = c.cellvalues

"""
    Ferrite.celldofs(c::CellBuffer)

`Ferrite.jl`'s `celldofs` function is overloaded on the `CellBuffer` to return 
the current cell's degree of freedom indices.
"""
@inline Ferrite.celldofs(c::CellBuffer) = c.dofs

# Convenience access functions to CellBuffer for use inside element routines
"""
    Ferrite.getcoordinates(c::CellBuffer)

`Ferrite.jl`'s `getcoordinates` function is overloaded on the `CellBuffer` to return 
the current cell's nodal coordinates. 
"""
@inline Ferrite.getcoordinates(c::CellBuffer) = c.coords

"""
    FerriteAssembly.get_aeold(c::CellBuffer)

Get the old element dof-values for the current cell 
(Filled by `NaN`s unless `aold` is passed to `doassemble!`)
"""
@inline get_aeold(c::CellBuffer) = c.ae_old

"""
    Ferrite.cellid(c::CellBuffer)

Get the current cell id/nr
"""
@inline Ferrite.cellid(c::CellBuffer) = c.cellid

"""
    Ferrite.dof_range(c::CellBuffer, name::Symbol)

Get the `dofrange::UnitRange{Int}` for the dofs pertaining to the field: `name`.
Same output as dof_range(dh::DofHandler, name), but fully type-stable. 
"""
@inline Ferrite.dof_range(c::CellBuffer, name::Symbol) = c.dofrange[name]

Ferrite.getfieldnames(c::CellBuffer) = keys(c.dofrange)

"""
    get_state_old(c::CellBuffer)

Get the state variables for the cell from the previous time step. 
"""
@inline get_state_old(c::CellBuffer) = c.old_state

"""
    get_time_increment(c::CellBuffer)

Get the time increment, `Δt`, that was passed to `doassemble`
"""
@inline get_time_increment(c::CellBuffer) = c.Δt

"""
    FerriteAssembly.get_user_data(c::CellBuffer)

Get the user specified `user_data` given to `CellBuffer`
"""
@inline get_user_data(c::CellBuffer) = c.user_data

"""
    FerriteAssembly.get_cache(c::CellBuffer)

Get the user-specified `cache` given to the `CellBuffer`
"""
@inline get_cache(c::CellBuffer) = c.cache

"""
    Ferrite.reinit!(c::CellBuffer, dh::AbstractDofHandler, cellnum::Int, anew, aold)

Reinitialize the `c::CellBuffer` for cell number `cellnum`.
The global degree of freedom vectors `anew` (current) and `aold` are used
to update the cell degree of freedom vectors in `c`.
If the global vectors are instead `::Nothing`, the corresponding cell values are set to `NaN`
The element stiffness, `c.Ke`, and residual, `c.re`, are also zeroed. 
"""
function Ferrite.reinit!(c::CellBuffer, dh::Ferrite.AbstractDofHandler, cellnum::Int, anew, aold, old_states)
    c.cellid = cellnum
    # `update_old_state!` could be made overloadable to make it error for certain state types 
    update_old_state!(::Any, ::Nothing) = nothing # Allow nothing to avoid having to pass old state if not needed. 
    update_old_state!(cb, _old_states) = (cb.old_state = fast_getindex(_old_states, cellnum))
    update_old_state!(c, old_states)
    celldofs!(c.dofs, dh, cellnum)
    getcoordinates!(c.coords, dh.grid, cellnum)
    reinit!(c.cellvalues, c.coords)
    _copydofs!(c.ae, anew, c.dofs)     # ae_new .= a_new[dofs]
    _copydofs!(c.ae_old, aold, c.dofs) # ae_old .= a_old[dofs]
    fill!(c.Ke, 0)
    fill!(c.re, 0)
    return nothing  # Ferrite's reinit! doesn't return 
end

# End of required functions for a custom CellBuffer