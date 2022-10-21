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
        cellvalues::Union{CellValues,Nothing}=nothing)

Creates a `Vector` of states for each cell in `dh`, where a `state=statefun(x::Vec)` and `x` is 
the coordinate in the grid. If `isnothing(cellvalues)`, then the `statefun` is called once for 
the average nodal coordinate of each cell. Otherwise, it is called once per quadrature point in 
`cellvalues` for that cell. 
"""
function create_states(dh::DofHandler, statefun::Function=Returns(nothing), cellvalues=nothing)
    return map(cell->create_state(statefun, cell, cellvalues), CellIterator(dh))
end

"""
    create_states(dh::MixedDofHandler, 
        statefun::Union{Function,Tuple}=Returns(nothing), 
        cellvalues::Union{CellValues,Tuple,Nothing}=nothing
        )

Creates a `Tuple` of the output of the following function for each fieldhandler in `dh`.
If `statefuns` and/or `cellvalues` are not tuples, the same value is used for each fieldhandler. 

    create_states(dh::MixedDofHandler, fh::FieldHandler, 
        statefun::Function, cellvalues::Union{Nothing,CellValues})

Returns a `Dict{Int}` with states for each cell in `fh`'s `cellset` and keys corresponding to the global `cellid`
If `isnothing(cellvalues)`, then the `statefun(x::Vec)` is called once for the average nodal coordinate `x` of each cell. 
Otherwise, it is called once for each quadrature point location `x` given by `cellvalues` for that cell. 
"""
function create_states(dh::MixedDofHandler, statefuns=Returns(nothing), cellvalues=nothing)
    num_fh = length(dh.fieldhandlers)
    statefuns_ = _maketuple(statefuns, num_fh)
    cellvalues_ = _maketuple(cellvalues, num_fh)
    return ntuple(i->create_states(dh, dh.fieldhandlers[i], statefuns_[i], cellvalues_[i]), num_fh)
end
function create_states(dh::MixedDofHandler, fh::FieldHandler, statefun::Function, cellvalues::Union{Nothing,CellValues})
    return Dict(cellid(cell)=>create_state(statefun, cell, cellvalues) for cell in CellIterator(dh, collect(fh.cellset)))
end