"""
    create_cell_state(material, cellvalues, x, ae)

Defaults to returning `nothing`.

Overload this function to create the state which should be passed 
into the `element_routine!`/`element_residual!` for the given 
material and cellvalues. 
As for the element routines, `ae`, is filled with `NaN` unless the
global degree of freedom vector is given to the 
[`create_states`](@ref) function. 
"""
create_cell_state(args...) = nothing

"""
    _create_cell_state(cell, material, cellvalues, a, ae)

Internal function to extract the relevant quantities from the 
`CellIterator`, `cell`: Pass these into the 
`create_cell_state` function that the user should specify. 
"""
function _create_cell_state(cell, material, cellvalues, a, ae)
    x = getcoordinates(cell)
    isnothing(cellvalues) || reinit!(cellvalues, x)
    _copydofs!(ae, a, celldofs(cell))
    return create_cell_state(material, cellvalues, x, ae)
end


"""
    create_states(dh::DofHandler, material=nothing, cellvalues=nothing, a=nothing)

Create the state variables for the given `dh`, `material`, `cellvalues`, and global 
dof vector `a`. `material` and `cellvalues` can be a Dict{String}, in which case they 
should have the same keys, corresponding to cellsets in the grid. 
This follows the same structure as for creating a `CellBuffer`.
The function `create_cell_state` is called for each cell. 
"""
function create_states(dh::DofHandler, material=nothing, cellvalues=nothing, a=nothing)
    ae = zeros(length(celldofs(dh, 1)))
    return map(cell->_create_cell_state(cell, material, cellvalues, a, ae), CellIterator(dh))
end

"""
    create_states(dh::MixedDofHandler, material=nothing, cellvalues=nothing, a=nothing)

Create the state variables for the given `dh`, `material`, `cellvalues`, and global dof vector `a`. 

`material` and `cellvalues` follow the same input structure as to `setup_cellbuffer`. 
I.e., they can be tuples with the same length as `fh.fieldhandlers`. 
Note that if `cellvalues::Tuple` is desired to be passed into the `create_cell_state`
function, it must be wrapped in an outer ntuple with the same length as `fh.fieldhandlers`.
(But using a `NamedTuple` is recommended to avoid this problem). 

For multiple materials, `material::Dict{String}`, and `cellvalues` can optionally be that as well, 
but then they must have the same keys. This follows the same structure as for creating a `CellBuffer`.
The function `create_cell_state` is called for each cell. 

The order of each call to `create_cell_state` is fixed for a given input, i.e. the loops are only 
over sorted `Set`s and `Dict` keys. (Useful if using `StableRNGs.jl` for random variations to obtain
reproducible results)
"""
function create_states(dh::MixedDofHandler, materials=nothing, cellvalues=nothing, a=nothing, cellset=nothing)
    return create_states_for_set(dh::MixedDofHandler, materials, cellvalues, a, cellset)
end

# Special dispatches for multiple materials in each DofHandler/FieldHandler
# Returns Dict{String,Dict{Int}} with cellsets as "outer" keys and cellid as "inner" keys
function create_states(dh::DofHandler, materials::Dict{String}, cellvalues=nothing, a=nothing)
    return _dict_create_states(dh, materials, cellvalues, a)
end
function create_states(dh::MixedDofHandler, materials::Dict{String}, cellvalues=nothing, a=nothing)
    return _dict_create_states(dh, materials, cellvalues, a)
end
function _dict_create_states(dh, materials, cellvalues, a)
    setnames = keys(materials)
    _cellvalues = _makedict(cellvalues, setnames)
    return Dict(
        setname=>create_states_for_set(dh, materials[setname], _cellvalues[setname], a, getcellset(dh,setname))
        for setname in sort(collect(setnames)))
end

# Functions that return a Dict{Int,State}
function create_states_for_set(dh::DofHandler, material, cellvalues, a, cellset)
    ae = zeros(length(celldofs(dh, first(cellset))))    # Here also ok to use (dh, 1), could use Ferrite.ndofs_per_cell
    return Dict(cellid(cell)=>_create_cell_state(cell, material, cellvalues, a, ae) for cell in CellIterator(dh, sort(collect(cellset))))
end
function create_states_for_set(dh::MixedDofHandler, materials, cellvalues, a, cellset)
    num_fh = length(dh.fieldhandlers)
    materials_ = _maketuple(materials, num_fh)
    cellvalues_ = _maketuple(cellvalues, num_fh)
    return ntuple(i->create_states_fh(dh, dh.fieldhandlers[i], materials_[i], cellvalues_[i], a, cellset), num_fh)
end

function create_states_fh(dh::MixedDofHandler, fh::FieldHandler, material, cellvalues, a, cellset)
    set = sort(collect(intersect_nothing(fh.cellset, cellset)))
    ae = zeros(length(celldofs(dh, first(set))))    # Could use Ferrite.ndofs_per_cell, but not documented
    return Dict(cellid(cell)=>_create_cell_state(cell, material, cellvalues, a, ae) for cell in CellIterator(dh, set))
end
