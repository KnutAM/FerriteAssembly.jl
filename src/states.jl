"""
    create_cell_state(material, cellvalues, x, ae, dofrange)

Defaults to returning `nothing`.

Overload this function to create the state which should be passed into the 
`element_routine!`/`element_residual!` for the given `material` and `cellvalues`. 
`x` is the cell's coordinates, `ae` the element degree of freedom values, and 
`dofrange::NamedTuple` containing the local dof range for each field. 
As for the element routines, `ae`, is filled with `NaN` unless the global degree 
of freedom vector is given to the [`setup_assembly`](@ref) function.
"""
create_cell_state(args...) = nothing
create_cell_state(material, cellvalues, x, ae, dh_fh) = create_cell_state(material, cellvalues, x, ae) # Backwards comp

"""
    _create_cell_state(cell, material, cellvalues, a, ae, dofrange, cellnr)

Internal function to reinit and extract the relevant quantities from the 
`cell::CellCache`, reinit cellvalues, update `ae` from `a`, and 
pass these into the `create_cell_state` function that the user should specify. 
"""
function _create_cell_state(cell, material, cellvalues, a, ae, dofrange, cellnr)
    reinit!(cell, cellnr)
    x = getcoordinates(cell)
    reinit!(cellvalues, x)
    _copydofs!(ae, a, celldofs(cell))
    return create_cell_state(material, cellvalues, x, ae, dofrange)
end

"""
    create_states(sdh::SubDofHandler, material, cellvalues, a, cellset, dofrange)

Create a `Dict` of states for the cells in `cellset`, where the user should 
define the [`create_cell_state`](@ref) function for their `material` (and corresponding `cellvalues`)
`dofrange::NamedTuple` is passed onto `create_cell_state` and contains the local dof ranges for each field. 
"""
function create_states(sdh::SubDofHandler, material, cellvalues, a, cellset, dofrange)
    ae = zeros(ndofs_per_cell(sdh))
    cell = CellCache(sdh.dh)
    return Dict(cellnr => _create_cell_state(cell, material, cellvalues, a, ae, dofrange, cellnr) for cellnr in cellset)
end

"""
    update_states!(old_states, new_states)

Update `old_states` with the values from `new_states`. This is typically done after a converged time increment.

This method tries to avoid allocating new values where possible. 
Currently, if [`create_cell_state`](@ref) returns `T` or `Vector{T}` where `isbitstype(T)`, this works.

If needed/wanted, it should be relatively easy to provide an interface to make it possible to have allocation free 
for custom cell states.
""" 
function update_states!(old_states::T, new_states::T) where T<:Union{Dict{String,<:Dict{Int}}, Dict{Int,<:Vector}}
    for (key, new_s) in new_states
        update_states!(old_states[key], new_s)
    end
end
update_states!(::T, ::T) where T<:Union{Vector{Nothing},Dict{Int,Nothing}} = nothing 

@inline function update_states!(old_states::T, new_states::T) where T<:Union{Vector{ET},Dict{Int,ET}} where ET
    copy_states!(Val(isbitstype(ET)), old_states, new_states)
end
@inline copy_states!(::Val{true},  old_states, new_states) = copy!(old_states,new_states)
@inline copy_states!(::Val{false}, old_states, new_states) = copy!(old_states,deepcopy(new_states))

#= # Something like this with replacing old_states[key] = deepcopy(new_s) with 
# copy_states!(old_states, key, new_s) to dispatch on new_s could work, 
# but otherwise directy copy!(dst, deepcopy(src)) seems better
function copy_states!(::Val{false}, old_states, new_states)
    for (key, new_s) in pairs(new_states)
        old_states[key] = deepcopy(new_s)
    end
end
=#

#= # Not used anymore (required if states are stored as vectors at top level instad of Dict{Int})
function update_states!(old_states::T, new_states::T) where T<:Vector{<:Vector}
    for (old_s, new_s) in zip(old_states, new_states)
        update_states!(old_s, new_s)
    end
end
=#