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
