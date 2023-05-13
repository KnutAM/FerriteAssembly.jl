"""
    AssemblyDomain(name, dh, material, cellvalues; cellset, colors=nothing, user_data=nothing)

Create an `AssemblyDomain` that can be used when calling [`setup_assembly`](@ref) to assemble multiple domains.
`name` is used to access the corresponding `DomainBuffer` and state variables returned by `setup_assembly`. 
If not given, `cellset` is attempted to be inferred from the DofHandler, `dh`. 
"""
struct AssemblyDomain
    name::String
    sdh::SubDofHandler
    material::Any
    cellvalues::Union{CellValues,NamedTuple}
    cellset::Union{AbstractVector{Int},AbstractSet{Int}} # Includes UnitRange
    colors::Any
    user_data::Any
end
function AssemblyDomain(name, sdh, material, cellvalues; cellset=getcellset(sdh), colors=nothing, user_data=nothing)
    return AssemblyDomain(name, sdh, material, cellvalues, cellset, colors, user_data)
end
function AssemblyDomain(name, dh::DofHandler, args...; kwargs...)
    return AssemblyDomain(name, SubDofHandler(dh), args...; kwargs...)
end

abstract type AbstractDomainBuffer end

struct DomainBuffer{SDH<:SubDofHandler,CB<:AbstractCellBuffer} <: AbstractDomainBuffer
    # TA=true if threaded, false if sequential
    sdh::SDH                # SubDofHandler
    cellset::Vector{Int}    # 
    cellbuffer::CB          # AbstractCellBuffer
end
function DomainBuffer(sdh::SubDofHandler, cellset, cellbuffer::AbstractCellBuffer)
    cellset_ = sort!(collect(intersect(getcellset(sdh), cellset))) # Saves a copy that is also sorted, making sure it is always Vector{Int}
    return DomainBuffer(sdh, cellset_, cellbuffer)
end

struct ThreadedDomainBuffer{SDH<:SubDofHandler, CB<:AbstractCellBuffer} <: AbstractDomainBuffer
    sdh::SDH 
    cellset::Vector{Int}
    colors::Vector{Vector{Int}}
    cellbuffers::TaskLocals{CB,CB}
end
# Threaded
function ThreadedDomainBuffer(sdh::SubDofHandler, cellset, colors::Vector, cellbuffer::AbstractCellBuffer)
    cellset_intersect = sort!(collect(intersect(cellset, getcellset(sdh))))
    colors_intersect = map(sort! ∘ collect ∘ Base.Fix1(intersect, cellset_intersect), colors)
    cellbuffers = TaskLocals(cellbuffer)
    return ThreadedDomainBuffer(sdh, cellset_intersect, colors_intersect, cellbuffers)
end

Ferrite.getcellset(b::Union{DomainBuffer, ThreadedDomainBuffer}) = b.cellset

# Type-unstable switch
setup_domainbuffer(threaded::Bool, args...; kwargs...) = setup_domainbuffer(Val(threaded), args...; kwargs...)

# Sequential
setup_domainbuffer(::Val{false}, sdh, cellset, cellbuffer, colors) = DomainBuffer(sdh, cellset, cellbuffer)

# Threaded
function setup_domainbuffer(::Val{true}, sdh, cellset, cellbuffer, colors::Nothing)
    _colors = create_coloring(_getgrid(sdh), cellset)
    return ThreadedDomainBuffer(sdh, cellset, _colors, cellbuffer)
end
function setup_domainbuffer(::Val{true}, sdh, cellset, cellbuffer, colors)
    return ThreadedDomainBuffer(sdh, cellset, colors, cellbuffer)
end

"""
    setup_assembly(dh, material, cellvalues; kwargs...)

Setup assembly for a single domain (i.e. the same material and interpolations everywhere).
Returns the `buffer`, `old_states`, and `new_states` to be used in [`doassemble!`](@ref).
The state variables are created by calling the user-defined [`create_cell_state`](@ref) function. 
Available keyword arguments
- `a=nothing`: Give the global dof vector to pass element dofs, `ae`, to `create_cell_state` (`NaN`-values otherwise)
- `threaded=Val(false)`: Set to `Val(true)` to setup threaded assembly. More elaborate settings may be added in the future. 
- `autodiffbuffer=Val(false)`: Set to `true` or `Val(true)` (for type stable construction) to use `AutoDiffCellBuffer`
  instead of `CellBuffer`.
- `user_data=nothing`: The `user_data` is passed to each `AbstractCellBuffer` by reference (when threaded)
- `colors=nothing`: Give colors for the grid from `Ferrite.create_coloring` to setup threaded colored assembly. 
  If `nothing`, Ferrite's default coloring algorithm is used.
- `cellset="all cells in dh"`: Which cells to assemble. In most cases, it is better to setup `AssemblyDomain`s 
  to assemble different domains. But this option can be used to assemble only a subset of the grid.
"""
setup_assembly(dh::DofHandler, args...; kwargs...) = setup_assembly(SubDofHandler(dh), args...; kwargs...)
function setup_assembly(sdh::SubDofHandler, material, cellvalues; a=nothing, user_data=nothing, 
        autodiffbuffer=Val(false), threading=Val(false), colors=nothing, cellset=getcellset(sdh)
        )
    dofrange = create_dofrange(sdh)
    new_states = create_states(sdh, material, cellvalues, a, cellset, dofrange)
    old_states = create_states(sdh, material, cellvalues, a, cellset, dofrange)

    cell_state = last(first(old_states))
    cellbuffer = setup_cellbuffer(autodiffbuffer, sdh, cellvalues, material, cell_state, dofrange, user_data)

    domainbuffer = setup_domainbuffer(threading, sdh, cellset, cellbuffer, colors)

    return domainbuffer, new_states, old_states
end

"""
    setup_assembly(domains::Vector{<:AssemblyDomain}; kwargs...)

Setup assembly for each [`AssemblyDomain`](@ref) in `domains`. 
Returns the `Dict`s `buffers`, `old_states`, and `new_states` to be used in [`doassemble!`](@ref).
The state variables are created by calling the user-defined [`create_cell_state`](@ref) function. 
Available keyword arguments
- `a=nothing`: Give the global dof vector to pass element dofs, `ae`, to `create_cell_state` (`NaN`-values otherwise)
- `autodiffbuffer=Val(false)`: Set to `true` or `Val(true)` (for type stable construction) to use `AutoDiffCellBuffer`
  instead of `CellBuffer`.
- `threading=Val(false)`: Set to `Val(true)` to setup threaded assembly. More elaborate settings may be added in the future. 
"""
function setup_assembly(domains::Vector{<:AssemblyDomain}; a=nothing, autodiffbuffer=Val(false), threading=Val(false))
    buffers = Dict{String,Any}()
    new_states = Dict{String,Dict{Int}}()
    old_states = Dict{String,Dict{Int}}()
    num_cells_grid = getncells(_getgrid(domains[1].sdh))
    added_cells = sizehint!(Set{Int}(), num_cells_grid)
    local buffer # Gives access after loop
    for d in domains
        n = d.name
        buffer, new_states[n], old_states[n] =
            setup_assembly(d.sdh, d.material, d.cellvalues; 
                a=a, threading=threading, autodiffbuffer=autodiffbuffer, 
                cellset=d.cellset, colors=d.colors, user_data=d.user_data)
        buffers[n] = buffer
        
        # Checks if some cells have been added before, and warn if that is the case.
        old_num = length(added_cells)
        num_add = length(getcellset(buffer))
        union!(added_cells, getcellset(buffer))
        length(added_cells) < (old_num + num_add) && @warn("The domain $n has cells that overlap with previous domains")
    end
    # Check that all cells have been added and warn otherwise. 
    num_cells_added = length(added_cells)
    num_cells_grid != num_cells_added && @warn("There are $num_cells_grid cells in the grid, but only $num_cells_added have been added")
    # Make the returned buffer dict have the type of buffer as a type parameter (for dispatch)
    BT = isa(buffer, DomainBuffer) ? DomainBuffer : ThreadedDomainBuffer
    buffers_typed = Dict{String,BT}(key=>val for (key,val) in buffers) # 
    return buffers_typed, new_states, old_states
end