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
function create_states(dh::DofHandler, statefun::Function=Returns(nothing), cellvalues=nothing)
    return map(cell->create_state(statefun, cell, cellvalues), CellIterator(dh))
end
# Not documented as only used internally to avoid clottering documentation with implementation details. 
function create_states(dh::DofHandler, statefun::Function, cellvalues, cellset)
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