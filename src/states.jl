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
create_cell_state(args...) = nothing

"""
    _create_cell_state(cell, material, cellvalues, a, ae, dofrange, cellnr)

Internal function to reinit and extract the relevant quantities from the 
`cell::CellCache`, reinit cellvalues, update `ae` from `a`, and 
pass these into the `create_cell_state` function that the user should specify. 
"""
function _create_cell_state(coords, dofs, material, cellvalues, a, ae, dofrange, sdh, cellnr)
    getcoordinates!(coords, _getgrid(sdh), cellnr)
    reinit!(cellvalues, coords)
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

"""
    update_states!(old_states, states)

In most cases, this 2-argument function is not required, and the 
entire domain buffer can be passed instead, see
[`update_states!`](@ref update_states!(::FerriteAssembly.DomainBuffers)).
This 2-argument function will then be called for the stored state variables. 
""" 
function update_states!(old_states::T, states::T) where T<:Dict{Int,<:Vector}
    for (key, new_s) in states
        update_states!(old_states[key], new_s)
    end
end
function update_states!(::T, ::T) where T<:Union{Vector{Nothing},Dict{Int,Nothing},Dict{Int,Vector{Nothing}}}
    return nothing 
end
@inline function update_states!(old_states::T, states::T) where T<:Union{Vector{ET},Dict{Int,ET}} where ET
    copy_states!(Val(isbitstype(ET)), old_states, states)
end
@inline copy_states!(::Val{true},  old_states, states) = copy!(old_states,states)
@inline copy_states!(::Val{false}, old_states, states) = copy!(old_states,deepcopy(states))
