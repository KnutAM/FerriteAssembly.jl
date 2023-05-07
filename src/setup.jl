"""
    AssemblyDomain(name, dh, material, cellvalues; cellset, user_data=nothing, cache=nothing)

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
    user_data::Any
    cache::Any
end
function AssemblyDomain(name, sdh, material, cellvalues; cellset=getcellset(sdh), user_data=nothing, cache=nothing)
    return AssemblyDomain(name, sdh, material, cellvalues, cellset, user_data, cache)
end
function AssemblyDomain(name, dh::DofHandler, args...; kwargs...)
    return AssemblyDomain(name, SubDofHandler(dh), args...; kwargs...)
end

struct DomainBuffer{TA,SDH<:SubDofHandler,CS,CB,SC}
    # TA=true if threaded, false if sequential
    sdh::SDH        # SubDofHandler
    cellset::CS     # Vector{Vector{Int}} or Vector{Int} (threaded vs sequential)
    cellbuffer::CB  # Vector{AbstractCellBuffer} or AbstractCellBuffer (threaded vs sequential)
    scaling::SC     # Note: All `DomainBuffer`s share the same `scaling` and `scalings`
    scalings::Vector{SC} # For threaded assembly (empty for sequential)
    function DomainBuffer{TA}(sdh::SubDofHandler, cellset, cellbuffer, scaling::SC, scalings::Vector{SC}) where {TA, SC}
        @assert isa(TA, Bool)
        return new{TA, typeof(sdh), typeof(cellset), typeof(cellbuffer), SC}(sdh, cellset, cellbuffer, scaling, scalings)
    end
end
# Sequential
function DomainBuffer(#=colors=#::Nothing, sdh::SubDofHandler, cellset, cellbuffer::AbstractCellBuffer, scaling::SC) where SC
    cellset_ = sort!(collect(cellset)) # Saves a copy that is also sorted, making sure it is always Vector{Int}
    return DomainBuffer{false}(sdh, cellset_, cellbuffer, scaling, SC[])
end
# Threaded
function DomainBuffer(colors::Vector, sdh::SubDofHandler, cellset, cellbuffer::AbstractCellBuffer, scaling)
    cellsets = map(sort! ∘ collect ∘ Base.Fix1(intersect, cellset), colors)
    cellbuffers = [copy_for_threading(cellbuffer) for _ in 1:Threads.nthreads()]
    # scaling will be Tuple{T, Vector{T}} if we have multiple domains, in which case distribute_scalings just returns scaling
    # For single domains, scaling will just be T where T is the type of scaling, and distribute_scalings returns T, Vector{T}
    global_scaling, perthread_scalings = distribute_scalings(scaling, colors)
    return DomainBuffer{true}(sdh, cellsets, cellbuffers, global_scaling, perthread_scalings)
end

# iscolored(::DomainBuffer{TA}) where TA = TA
iscolored(::Dict{String,<:DomainBuffer{TA}}) where TA = TA

"""
    setup_assembly(dh, material, cellvalues; kwargs...)

Setup assembly for a single domain (i.e. the same material and interpolations everywhere).
Returns the `buffer`, `old_states`, and `new_states` to be used in [`doassemble!`](@ref).
The state variables are created by calling the user-defined [`create_cell_state`](@ref) function. 
Available keyword arguments
- `a=nothing`: Give the global dof vector to pass element dofs, `ae`, to `create_cell_state` (`NaN`-values otherwise)
- `colors=nothing`: Give colors for the grid from `Ferrite.create_coloring` to setup threaded assembly. 
  If `nothing`, the assembly is sequential.
- `autodiffbuffer=Val(false)`: Set to `true` or `Val(true)` (for type stable construction) to use `AutoDiffCellBuffer`
  instead of `CellBuffer`.
- `user_data=nothing`: The `user_data` is passed to each `AbstractCellBuffer` by reference (when threaded)
- `cache=nothing`: The `cache` is passed to each `AbstractCellBuffer`, and deepcopied if threaded. 
- `scaling=nothing`: The scaling to be calculated, see e.g. [`ElementResidualScaling`](@ref)
"""
setup_assembly(dh::DofHandler, args...; kwargs...) = setup_assembly(SubDofHandler(dh), args...; kwargs...)
function setup_assembly(sdh::SubDofHandler, material, cellvalues; cellset=getcellset(sdh), a=nothing,
        colors=nothing, autodiffbuffer=Val(false), user_data=nothing, cache=nothing, scaling=nothing
        )
    dofrange = create_dofrange(sdh)
    new_states = create_states(sdh, material, cellvalues, a, cellset, dofrange)
    old_states = create_states(sdh, material, cellvalues, a, cellset, dofrange)

    cell_state = last(first(old_states))
    cellbuffer = setup_cellbuffer(autodiffbuffer, sdh, cellvalues, material, cell_state, dofrange, user_data, cache)

    domainbuffer = DomainBuffer(colors, sdh, cellset, cellbuffer, scaling)

    return domainbuffer, old_states, new_states
end

"""
    setup_assembly(domains::Vector{<:AssemblyDomain}; kwargs...)

Setup assembly for each [`AssemblyDomain`](@ref) in `domains`. 
Returns the `Dict`s `buffers`, `old_states`, and `new_states` to be used in [`doassemble!`](@ref).
The state variables are created by calling the user-defined [`create_cell_state`](@ref) function. 
Available keyword arguments
- `a=nothing`: Give the global dof vector to pass element dofs, `ae`, to `create_cell_state` (`NaN`-values otherwise)
- `colors=nothing`: Give colors for the grid from `Ferrite.create_coloring` to setup threaded assembly. 
  If `nothing`, the assembly is sequential.
- `autodiffbuffer=Val(false)`: Set to `true` or `Val(true)` (for type stable construction) to use `AutoDiffCellBuffer`
  instead of `CellBuffer`.
- `scaling=nothing`: The scaling to be calculated, see e.g. [`ElementResidualScaling`](@ref)
"""
function setup_assembly(domains::Vector{<:AssemblyDomain}; colors=nothing, scaling=nothing, kwargs...)
    is_threaded = !isnothing(colors)
    buffers = Dict{String,DomainBuffer{is_threaded}}()
    new_states = Dict{String,Dict{Int}}()
    old_states = Dict{String,Dict{Int}}()
    scalings = distribute_scalings(scaling, colors)
    num_cells = 0
    all_cells = Set{Int}()
    for d in domains
        n = d.name
        buffers[n], old_states[n], new_states[n] = 
            setup_assembly(d.sdh, d.material, d.cellvalues; 
                cellset=d.cellset, user_data=d.user_data, cache=d.cache,
                colors=colors, scaling=scalings, kwargs...
                )
        num_cells += length(d.cellset)
        all_cells = union!(all_cells, d.cellset)
        length(all_cells) != num_cells && @warn("The domain $n has cells that overlap with previous domains")
    end
    num_cells_grid = getncells(domains[1].sdh.dh.grid)
    num_cells = length(all_cells)
    num_cells != num_cells_grid && @warn("There are $num_cells_grid cells in the grid, but only $num_cells have been added")
    return buffers, old_states, new_states
end

function distribute_scalings(scaling, ::Vector) # Distribute to each thread
    return (scaling, [deepcopy(scaling) for _ in 1:Threads.nthreads()]) 
end
distribute_scalings(scaling, ::Nothing) = scaling   # Sequential
# Already distributed:
distribute_scalings(scalings::Tuple{T,Vector{T}}, ::Vector) where T = scalings
