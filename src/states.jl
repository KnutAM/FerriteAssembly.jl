# Minimal interface for a vector, storage format will probably be updated later. 
mutable struct StateVector{SV}
    vals::Dict{Int, SV}
end
Base.getindex(s::StateVector, cellnum::Int) = s.vals[cellnum]
Base.setindex!(s::StateVector, v, cellnum::Int) = setindex!(s.vals, v, cellnum)
Base.:(==)(a::StateVector, b::StateVector) = (a.vals == b.vals)

struct StateVariables{SV}
    old::StateVector{SV} # Rule: Referenced during assembly, not changed (ever)
    new::StateVector{SV} # Rule: Updated during assembly, not referenced (before updated)
end
StateVariables(old::Dict, new::Dict) = StateVariables(StateVector(old), StateVector(new))

function update_states!(sv::StateVariables)
    tmp = sv.old.vals
    sv.old.vals = sv.new.vals
    sv.new.vals = tmp
end

"""
    copy_state(state)

Return a copy of `state` such that mutating the returned value does not affect `state`.

Used by [`set_new_to_old_states!`](@ref) to copy the old state into the new state, for
cell states that are not `AbstractArray`s (`AbstractArray` states are copied with
`copyto!` instead, which avoids allocations).

Defaults to `deepcopy(state)`. Overload this function for a custom cell state type
to provide a more efficient copy (e.g. by reusing already allocated memory), instead
of relying on the `deepcopy` fallback.
"""
copy_state(state) = deepcopy(state)

function set_new_to_old_states!(sv::StateVariables)
    for key in keys(sv.old.vals)
        old_val = sv.old.vals[key]
        if old_val isa AbstractArray
            copyto!(sv.new.vals[key], old_val)
        else
            sv.new.vals[key] = copy_state(old_val)
        end
    end
end

# Experimental, basically copy!, but use separate name for clarity
function replace_states!(dst::StateVariables, src::StateVariables)
    dst.old.vals = src.old.vals 
    dst.new.vals = src.new.vals
    return dst
end

"""
    create_cell_state(material, cellvalues, x, ae, dofrange)

Defaults to returning `nothing`.

Overload this function to create the state which should be passed into the 
`element_routine!`/`element_residual!` for the given `material` and `cellvalues`. 
`x` is the cell's coordinates, `ae` the element degree of freedom values, and 
`dofrange::NamedTuple` containing the local dof range for each field. 
As for the element routines, `ae`, is filled with `NaN` unless the global degree 
of freedom vector is given to the [`setup_domainbuffer`](@ref) function.
"""
create_cell_state(material, cv, args...) = [nothing for _ in 1:_getnquadpoints(cv)]

_getnquadpoints(fe_v::Ferrite.AbstractValues) = getnquadpoints(fe_v)
_getnquadpoints(nt::NamedTuple) = getnquadpoints(first(nt))

"""
    _create_cell_state(cell, material, cellvalues, a, ae, dofrange, cellnr)

Internal function to reinit and extract the relevant quantities from the 
`cell::CellCache`, reinit cellvalues, update `ae` from `a`, and 
pass these into the `create_cell_state` function that the user should specify. 
"""
function _create_cell_state(coords, dofs, material, cellvalues, a, ae, dofrange, sdh, cellnr)
    getcoordinates!(coords, _getgrid(sdh), cellnr)
    reinit!(cellvalues, getcells(_getgrid(sdh), cellnr), coords)
    celldofs!(dofs, sdh.dh, cellnr)
    _copydofs!(ae, a, dofs)
    return create_cell_state(material, cellvalues, coords, ae, dofrange)
end

"""
    create_states(sdh::SubDofHandler, material, cellvalues, a, cellset, dofrange)

Create a `Dict` of states for the cells in `cellset`, where the user should 
define the [`create_cell_state`](@ref) function for their `material` (and corresponding `cellvalues`)
`dofrange::NamedTuple` is passed onto `create_cell_state` and contains the local dof ranges for each field. 
"""
function create_states(sdh::SubDofHandler, material, cellvalues, a, cellset, dofrange)
    ae = zeros(ndofs_per_cell(sdh))
    coords = getcoordinates(_getgrid(sdh), first(cellset))
    dofs = zeros(Int, ndofs_per_cell(sdh))
    return Dict(cellnr => _create_cell_state(coords, dofs, material, cellvalues, a, ae, dofrange, sdh, cellnr) for cellnr in cellset)
end
