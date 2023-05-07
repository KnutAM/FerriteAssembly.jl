struct AssemblyDomain
    name::String
    sdh::SubDofHandler
    material::Any
    cellvalues::Union{CellValues,NamedTuple}
    cellset::OrderedSet{Int}
    user_data::Any
    cache::Any
end
function AssemblyDomain(name, sdh, material, cellvalues; cellset=getcellset(sdh), user_data=nothing, cache=nothing)
    return AssemblyDomain(name, sdh, material, cellvalues, cellset, user_data, cache)
end

struct DomainBuffer{TA,SDH,CS,CB,SC}
    # TA=true if threaded, false if sequential
    sdh::SDH        # SubDofHandler
    cellset::CS     # Vector{OrderedSet{Int}} or OrderedSet{Int} (threaded vs sequential)
    cellbuffer::CB  # Vector{AbstractCellBuffer} or AbstractCellBuffer (threaded vs sequential)
    scaling::SC     # Note: All `DomainBuffer`s share the same `scaling` and `scalings`
    scalings::Vector{SC} # For threaded assembly
    function DomainBuffer(sdh::SubDofHandler, cellset::OrderedSet, cellbuffer::AbstractCellBuffer, scaling, ::Nothing)
        return new{false, typeof(sdh), typeof(cellset), typeof(cellbuffer), typeof(scaling)}(sdh, cellset, cellbuffer, scaling, typeof(scaling)[])
    end
    function DomainBuffer(sdh::SubDofHandler, cellsets::Vector, cellbuffers::Vector, scalings::Tuple{SC, Vector{SC}}) where SC
        return new{true, typeof(sdh), typeof(cellsets), typeof(cellbuffers), SC}(sdh, cellsets, cellbuffers, scalings[1], scalings[2])
    end
end
function DomainBuffer(sdh::SubDofHandler, cellset::OrderedSet, cellbuffer::AbstractCellBuffer, scaling, colors::Vector)
    cellsets = map(Base.Fix1(intersect, cellset), colors)
    cellbuffers = [copy_for_threading(cellbuffer) for _ in 1:Threads.nthreads()]
    scalings = distribute_scalings(scaling, colors)
    return DomainBuffer(sdh, cellsets, cellbuffers, scalings)
end
iscolored(::Dict{String,<:DomainBuffer{true}}) = true
iscolored(::Dict{String,<:DomainBuffer{false}}) = false

"""
    setup_assembly(dh, material, cv_qr_quadorder; kwargs...)

Setup assembly for a single domain (i.e. the same material and interpolations everywhere).
Returns the `buffer`, `old_states`, and `new_states` to be used in [`doassemble!`](@ref).
The state variables are created by calling the user-defined [`create_cell_state`](@ref) function. 
"""
setup_assembly(dh::DofHandler, args...; kwargs...) = setup_assembly(SubDofHandler(dh), args...; kwargs...)
function setup_assembly(sdh::SubDofHandler, material, cv; cellset=getcellset(sdh), a=nothing,
        colors=nothing, autodiffbuffer=Val(false), user_data=nothing, cache=nothing, scaling=nothing
        )
    dofrange = create_dofrange(sdh)
    new_states = create_states(sdh, material, cv, a, cellset, dofrange)
    old_states = create_states(sdh, material, cv, a, cellset, dofrange)

    cell_state = last(first(old_states))
    cellbuffer = setup_cellbuffer(autodiffbuffer, sdh, cv, material, cell_state, dofrange, user_data, cache)

    domainbuffer = DomainBuffer(sdh, cellset, cellbuffer, scaling, colors)

    return domainbuffer, old_states, new_states
end

"""
    setup_assembly(domains::Vector{<:AssemblyDomain}; kwargs...)

Setup assembly for each domain in `domains`. 
Returns the `Dict`s `buffers`, `old_states`, and `new_states` to be used in [`doassemble!`](@ref).
The state variables are created by calling the user-defined [`create_cell_state`](@ref) function. 
"""
function setup_assembly(domains::Vector{<:AssemblyDomain}; colors=nothing, scaling=nothing, kwargs...)
    is_threaded = !isnothing(colors)
    buffers = Dict{String,DomainBuffer{is_threaded}}()
    new_states = Dict{String,OrderedDict{Int}}()
    old_states = Dict{String,OrderedDict{Int}}()
    scalings = distribute_scalings(scaling, colors)
    for d in domains
        n = d.name
        buffers[n], old_states[n], new_states[n] = 
            setup_assembly(d.sdh, d.material, d.cellvalues; 
                cellset=d.cellset, user_data=d.user_data, cache=d.cache,
                colors=colors, scaling=scalings, kwargs...
                )
    end
    return buffers, old_states, new_states
end

function distribute_scalings(scaling, ::Vector) # Distribute to each thread
    return (scaling, [deepcopy(scaling) for _ in 1:Threads.nthreads()]) 
end
distribute_scalings(scaling, ::Nothing) = scaling   # Sequential
# Already distributed:
distribute_scalings(scalings::Tuple{T,Vector{T}}, ::Vector) where T = scalings
