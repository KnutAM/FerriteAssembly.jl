struct AssemblyDomain
    name::String
    sdh::SubDofHandler
    material::Any
    cellvalues::Union{CellValues,NamedTuple}
    cellset::OrderedSet{Int}
    scaling::Any
    user_data::Any
    cache::Any
end
function AssemblyDomain(name, sdh, material, cellvalues; cellset=getcellset(sdh), scaling=nothing, user_data=nothing, cache=nothing)
    return AssemblyDomain(name, sdh, material, cellvalues, cellset, scaling, user_data, cache)
end

struct DomainBuffer{TA,SDH,CS,CB,SC}
    # TA=true if threaded, false if sequential
    sdh::SDH        # SubDofHandler
    cellset::CS     # Vector{OrderedSet{Int}} or OrderedSet{Int} (threaded vs sequential)
    cellbuffer::CB  # Vector{AbstractCellBuffer} or AbstractCellBuffer (threaded vs sequential)
    scaling::SC     # 
    scalings::Vector{SC} # For threaded assembly
    function DomainBuffer(sdh::SubDofHandler, cellset::OrderedSet, cellbuffer::AbstractCellBuffer, scaling, ::Nothing)
        return new{false, typeof(sdh), typeof(cellset), typeof(cellbuffer), typeof(scaling)}(sdh, cellset, cellbuffer, scaling, typeof(scaling)[])
    end
    function DomainBuffer(sdh::SubDofHandler, cellsets::Vector, cellbuffers::Vector, scaling::SC, scalings::Vector{SC}) where SC
        return new{true, typeof(sdh), typeof(cellsets), typeof(cellbuffers), typeof(scaling)}(sdh, cellsets, cellbuffers, scaling, scalings)
    end
end
function DomainBuffer(sdh::SubDofHandler, cellset::OrderedSet, cellbuffer::AbstractCellBuffer, scaling, colors::Vector)
    cellsets = map(Base.Fix1(intersect, cellset), colors)
    cellbuffers = [copy_for_threading(cellbuffer) for _ in 1:Threads.nthreads()]
    scalings = [deepcopy(scaling) for _ in 1:Threads.nthreads()]
    return DomainBuffer(sdh, cellsets, cellbuffers, scaling, scalings)
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
    setup_assembly(domains::Vector{<:AssemblyDomain}, a=nothing; kwargs...)

Setup assembly for each domain in `domains`. 
Returns the `Dict`s `buffers`, `old_states`, and `new_states` to be used in [`doassemble!`](@ref).
The state variables are created by calling the user-defined [`create_cell_state`](@ref) function. 
"""
function setup_assembly(domains::Vector{<:AssemblyDomain}, a=nothing; colors=nothing, kwargs...)
    is_threaded = !isnothing(colors)
    buffers = Dict{String,DomainBuffer{is_threaded}}()
    new_states = Dict{String,OrderedDict{Int}}()
    old_states = Dict{String,OrderedDict{Int}}()
    for d in domains
        n = d.name
        buffers[n], old_states[n], new_states[n] = 
            setup_assembly(d.sdh, d.material, d.cellvalues; cellset=d.cellset, a=a,
                colors=colors, user_data=d.user_data, cache=d.cache, scaling=d.scaling, kwargs...
                )
    end
    return buffers, old_states, new_states
end