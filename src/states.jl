# Minimal interface for a vector, storage format will probably be updated later. 
mutable struct StateVector{SV, VV <: AbstractVector{SV}}
    vals::VV
    inds::Vector{Int} # Can as an optimization be shared between all `StateVectors` (also across domains)
end
Base.getindex(s::StateVector, cellnum::Int) = s.vals[s.inds[cellnum]]
Base.setindex!(s::StateVector, v, cellnum::Int) = setindex!(s.vals, v, s.inds[cellnum])
Base.:(==)(a::StateVector, b::StateVector) = (a.vals == b.vals)
function Base.iterate(s::StateVector, i::Int = 1)
    i > length(s.vals) && return nothing
    return s.vals[i], i + 1
end

struct StateVariables{SV, VV}
    old::StateVector{SV, VV} # Rule: Referenced during assembly, not changed (ever)
    new::StateVector{SV, VV} # Rule: Updated during assembly, not referenced (before updated)
end
function StateVariables(inds::Vector{Int}, old::Dict{K, SV}, new::Dict{K, SV}) where {K, T, SV <: Vector{T}}
    # if eltype(old) isa Vector => ArrayOfVectorViews
    # else => Vector{eltype(old)}
    num = length(old)
    num_total = sum(length, values(old))
    old_data = Vector{T}(undef, num_total)
    new_data = Vector{T}(undef, num_total)
    indices = Vector{Int}(undef, num + 1)
    i = 1
    j = 1
    for key in sort(collect(keys(old)))
        indices[i] = j
        inds[key] = i
        for (oldval, newval) in zip(old[key], new[key])
            old_data[j] = oldval
            new_data[j] = newval
            j += 1
        end
        i += 1
    end
    indices[i] = j
    oldvals = ArrayOfVectorViews(indices, old_data, LinearIndices((num,)))
    newvals = ArrayOfVectorViews(indices, new_data, LinearIndices((num,)))
    StateVariables(StateVector(oldvals, inds), StateVector(newvals, inds))
end
function StateVariables(inds::Vector{Int}, old::Dict{K, SV}, new::Dict{K, SV}) where {K, SV}
    oldvals = Vector{SV}(undef, length(old))
    newvals = Vector{SV}(undef, length(new))
    i = 1
    for key in sort(collect(keys(old)))
        oldvals[i] = old[key]
        newvals[i] = new[key]
        inds[key] = i
        i += 1
    end
    StateVariables(StateVector(oldvals, inds), StateVector(newvals, inds))
end

function update_states!(sv::StateVariables)
    tmp = sv.old.vals
    sv.old.vals = sv.new.vals
    sv.new.vals = tmp
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
    grid = _getgrid(sdh)
    coords = getcoordinates(grid, first(cellset))
    dofs = zeros(Int, ndofs_per_cell(sdh))
    # Could make construction more efficient by doing this when creating the ArrayOfVectorViews
    old = Dict(cellnr => _create_cell_state(coords, dofs, material, cellvalues, a, ae, dofrange, sdh, cellnr) for cellnr in cellset)
    new = Dict(cellnr => _create_cell_state(coords, dofs, material, cellvalues, a, ae, dofrange, sdh, cellnr) for cellnr in cellset)
    inds = zeros(Int, getncells(grid)) # Could be moved out and shared between all domains...
    return StateVariables(inds, old, new)
end
