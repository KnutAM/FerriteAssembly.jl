"""
    create_cell_state(material, cellvalues, x, ae)

Overload this function to create the state which should be passed 
into the `element_routine!`/`element_residual!` for the given 
material and cellvalues. 
As for the element routines, `ae`, is filled with `NaN` if the 
global degree of freedom vector is given to the 
[`create_states`](@ref) function. 
If not defined, defaults to returning `nothing`
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
        for setname in setnames)
end

# Functions that return a Dict{Int,State}
function create_states_for_set(dh::DofHandler, material, cellvalues, a, cellset)
    ae = zeros(length(celldofs(dh, first(cellset))))    # Here also ok to use (dh, 1), could use Ferrite.ndofs_per_cell
    return Dict(cellid(cell)=>_create_cell_state(cell, material, cellvalues, a, ae) for cell in CellIterator(dh, collect(cellset)))
end
function create_states_for_set(dh::MixedDofHandler, materials, cellvalues, a, cellset)
    num_fh = length(dh.fieldhandlers)
    materials_ = _maketuple(materials, num_fh)
    cellvalues_ = _maketuple(cellvalues, num_fh)
    return ntuple(i->create_states_fh(dh, dh.fieldhandlers[i], materials_[i], cellvalues_[i], a, cellset), num_fh)
end

function create_states_fh(dh::MixedDofHandler, fh::FieldHandler, material, cellvalues, a, cellset)
    set = collect(intersect_nothing(fh.cellset, cellset))
    ae = zeros(length(celldofs(dh, first(set))))    # Could use Ferrite.ndofs_per_cell, but not documented
    return Dict(cellid(cell)=>_create_cell_state(cell, material, cellvalues, a, ae) for cell in CellIterator(dh, set))
end

"""
    create_states(dh::DofHandler, 
        statefun::Function=Returns(nothing), 
        cellvalues::Union{CellValues,Nothing}=nothing
        ) -> Vector

Creates a `Vector` of states for each cell in `dh`, where a `state=statefun(x::Vec)` and `x` is 
the coordinate in the grid. If `isnothing(cellvalues)`, then the `statefun` is called once for 
the average nodal coordinate of each cell. Otherwise, it is called once per quadrature point in 
`cellvalues` for that cell.

    create_states(dh::DofHandler, 
        statefun::Dict{String,Function}, 
        cellvalues::Union{Dict{String},CellValues,Nothing}}=nothing
        ) -> Dict{String,Dict{Int}}

    create_states(dh::MixedDofHandler, 
        statefun::Dict{String,Function}, 
        cellvalues::Union{Dict{String},CellValues,Nothing}}=nothing
        ) -> Dict{String,NTuple{N,Dict{Int}}}

If a `isa(statefun,Dict)`, then given a `Dict{String,Dict{Int}}` is returned,
one for each cellset in the grid corresponding to the keys of `statefun`. 
The elements are then `Dict`s that map the cellnr to the state for each cell 
in the cellset. If `cellvalues::Dict` is given, the keys should match 
`statefun`. Otherwise, the same cellvalues are assumed for each statefun.
This is typically used if multiple materials are used on the same grid.

    create_states(dh::MixedDofHandler, 
        statefun::Union{Function,Tuple}=Returns(nothing), 
        cellvalues::Union{CellValues,Tuple,Nothing}=nothing,
        cellset=nothing
        ) -> NTuple{N,Dict{Int}}

Creates a `Tuple` of the output of the following function for each fieldhandler in `dh`.
If `statefuns` and/or `cellvalues` are not tuples, the same value is used for each fieldhandler. 

    create_states(dh::MixedDofHandler, fh::FieldHandler, 
        statefun::Function, cellvalues::Union{Nothing,CellValues}, cellset
        ) -> Dict{Int}

Returns a `Dict{Int}` with states for each cell in `fh`'s `cellset` and keys corresponding to the global `cellid`
If `isnothing(cellvalues)`, then the `statefun(x::Vec)` is called once for the average nodal coordinate `x` of each cell. 
Otherwise, it is called once for each quadrature point location `x` given by `cellvalues` for that cell. 
If a `cellset!=nothing` is given, states are only created for the intersection of that `cellset` and the fieldhandlers `cellset` 
"""



#=
"""
    create_state(statefun, cell, cv::Nothing)
    create_state(statefun, cell, cv::CellValues)

Internal functions for calling `statefun(x)` for the cell center 
if cv::Nothing or for each quadrature point if cv::CellValues.
"""
function create_state(statefun, cell, ::Nothing)
    x = getcoordinates(cell)
    xc = sum(x)/length(x)
    return statefun(xc)
end
function create_state(statefun, cell, cv::CellValues)
    x = getcoordinates(cell)
    reinit!(cv, x)
    return [statefun(spatial_coordinate(cv, i, x)) for i in 1:getnquadpoints(cv)]
end

"""
    create_state(statefun, cell, cv, a, ae)

Internal function for calling `statefun(cv, coords, ae)` to create 
the state vector. Here, `coords` are the cell's nodal coordinates
"""
function create_state(statefun, cell, cv, a, ae)
    x = getcoordinates(cell)
    reinit!(cv, x)
    _copydofs!(ae, a, celldofs(cell))
    return statefun(cv, x, ae)
end

function create_states(dh::DofHandler, statefun::Function=Returns(nothing), cellvalues=nothing, a=nothing)
    ci = CellIterator(dh)
    ae = zeros(length(celldofs(ci)))
    return map(cell->create_state(statefun, cell, cellvalues, a, ae), CellIterator(dh))
end
# Not documented as only used internally to avoid clottering documentation with implementation details. 
function create_states_set(dh::DofHandler, statefun::Function, cellvalues, cellset)
    return Dict(cellid(cell)=>create_state(statefun, cell, cellvalues) for cell in CellIterator(dh, collect(cellset)))
end
create_states(dh::DofHandler, statefun::Dict, cellvalues=nothing) = _dict_create_states(dh, statefun, cellvalues)
create_states(dh::MixedDofHandler, statefun::Dict, cellvalues=nothing) = _dict_create_states(dh, statefun, cellvalues)
function _dict_create_states(dh::Ferrite.AbstractDofHandler, statefun::Dict, cellvalues)
    setnames = keys(statefun)
    _cellvalues = _makedict(cellvalues, setnames)
    return Dict(
        setname=>create_states(dh, statefun[setname], _cellvalues[setname], getcellset(dh,setname))
        for setname in setnames)
end

function create_states(dh::MixedDofHandler, statefuns=Returns(nothing), cellvalues=nothing, cellset=nothing)
    num_fh = length(dh.fieldhandlers)
    statefuns_ = _maketuple(statefuns, num_fh)
    cellvalues_ = _maketuple(cellvalues, num_fh)
    return ntuple(i->create_states(dh, dh.fieldhandlers[i], statefuns_[i], cellvalues_[i], cellset), num_fh)
end
function create_states(dh::MixedDofHandler, fh::FieldHandler, statefun::Function, cellvalues::Union{Nothing,CellValues}, cellset)
    set = collect(intersect_nothing(fh.cellset, cellset))
    return Dict(cellid(cell)=>create_state(statefun, cell, cellvalues) for cell in CellIterator(dh, set))
end
=#
