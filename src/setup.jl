"""
    DomainSpec(sdh::SubDofHandler, material, fe_values; set=getcellset(sdh), colors_or_chunks=nothing, user_data=nothing)
    DomainSpec(dh::DofHandler, material, fe_values; set=1:getncells(dh), colors=nothing, chunks=nothing, user_data=nothing)

Create a `DomainSpec` that can be used to set up a domain buffer. 

* `sdh`/`dh`: Give the `DofHandler` for the domain in question, or a `SubDofHandler` in case there are more than one in `DofHandler` 
  (See `Ferrite.jl's documentation`)
* `material`: Used for dispatch on the utilized `worker`'s function. 
* `fe_values`: `CellValues` or `FaceValues` depending on the type of domain. 
* `set`: The items in the domain, the element type determines the type of domain 
  - Cell domain: `Int`
  - Face domain: `FaceIndex`
* `colors::Vector{Vector{I}}`: used to avoid race conditions when multithreading. 
  For cell domains, `I=Int`, and for face domains, `I` can be either `Int` 
  (denoting cell numbers) or `FaceIndex` for actual faces. If `I=Int`, it will be 
  converted to `FaceIndex` internally. If `colors=nothing` and `chunks=nothing`, 
  `Ferrite.jl`'s default coloring algorithm is used.
* `chunks::Vector{Vector{Vector{I}}}`. During multithreading, each task works with 
  items in one `chunk::Vector{I}` at a time. Items in `chunks[k][i]` and `chunks[k][j]`
  should be independent (i.e. not share dofs). If given, this input takes precidence over 
  `colors`. For `chunks`, `I` must be `Int` for cell domains and `FaceIndex` for face domains. 
* `user_data`: Can be whatever the user wants to and is passed along by reference everywhere. 
  It is accessible from the `ItemBuffer` (e.g. `CellBuffer`) given to the `worker`'s function
  via the `get_user_data` function. However, since it is passed by reference, modifying values 
  during work, care must be taken to ensure thread safety. 
  To avoid allocations, caches can be created separately with `allocate_cell_cache` and 
  `allocate_face_cache`.
"""
struct DomainSpec{I}
    sdh::SubDofHandler
    material::Any
    fe_values::Any
    set::Union{AbstractVector{I},AbstractSet{I}}
    colors_or_chunks::Any
    user_data::Any
end
function DomainSpec(sdh::SubDofHandler, material, values; set=getcellset(sdh), colors=nothing, chunks=nothing, user_data=nothing)
    intersected_set = intersect_cellset_sort(set, getcellset(sdh))
    colors_or_chunks = chunks!==nothing ? chunks : colors
    return DomainSpec(sdh, material, values, intersected_set, colors_or_chunks, user_data)
end
function DomainSpec(dh::DofHandler, args...; kwargs...)
    length(dh.subdofhandlers) > 1 && error("A domain requires a single SubDofHandler")
    return DomainSpec(dh.subdofhandlers[1], args...; kwargs...)
end

"""
    setup_domainbuffers(domains::Dict{String,DomainSpec}; kwargs...)

Setup multiple domain buffers, one for each `DomainSpec` in `domains`.
See [`setup_domainbuffer`](@ref) for description of the keyword arguments.
"""
function setup_domainbuffers(domains::Dict{String,<:DomainSpec}; kwargs...)
    return Dict(name => setup_domainbuffer(domain; kwargs...) for (name, domain) in domains)
end

"""
    setup_domainbuffer(domain::DomainSpec; a=nothing, threading=false, autodiffbuffer=false)

Setup a domain buffer for a single grid domain, `domain`. 
* `a::Vector`: The global degree of freedom values are used to pass the 
  local element dof values to the [`create_cell_state`](@ref) function, 
  making it possible to create the initial state dependent on the initial 
  conditions for the field variables. 
* `threading`: Should a `ThreadedDomainBuffer` be created to work the grid multithreaded
  if supported by the used `worker`?
* `autodiffbuffer`: Should a custom itembuffer be used to speed up the automatic 
  differentiation (if supported by the itembuffer)
"""
function setup_domainbuffer(domain::DomainSpec; threading=Val(false), kwargs...)
    return _setup_domainbuffer(threading, domain; kwargs...)
end

create_states(domain::DomainSpec{Int}, a) = create_states(domain.sdh, domain.material, domain.fe_values, a, domain.set, create_dofrange(domain.sdh))
create_states(::DomainSpec{FaceIndex}, ::Any) = Dict{Int,Nothing}()

function setup_itembuffer(adb, domain::DomainSpec{FaceIndex}, args...)
    dofrange = create_dofrange(domain.sdh)
    return setup_facebuffer(adb, domain.sdh, domain.fe_values, domain.material, dofrange, domain.user_data)
end
function setup_itembuffer(adb, domain::DomainSpec{Int}, states)
    dofrange = create_dofrange(domain.sdh)
    return setup_cellbuffer(adb, domain.sdh, domain.fe_values, domain.material, first(values(states)), dofrange, domain.user_data)
end

function _setup_domainbuffer(threaded, domain; a=nothing, autodiffbuffer=Val(false))
    states = create_states(domain, a)
    old_states = create_states(domain, a)
    itembuffer = setup_itembuffer(autodiffbuffer, domain, states)
    return _setup_domainbuffer(threaded, domain.set, itembuffer, states, old_states, domain.sdh, domain.colors_or_chunks)
end

# Type-unstable switch
function _setup_domainbuffer(threaded::Bool, args...)
    return _setup_domainbuffer(Val(threaded), args...)
end
# Sequential
function _setup_domainbuffer(::Val{false}, set, itembuffer, states, old_states, sdh, args...)
    return DomainBuffer(set, itembuffer, states, old_states, sdh)
end
# Threaded
function _setup_domainbuffer(::Val{true}, set, itembuffer, states, old_states, sdh, args...)
    return ThreadedDomainBuffer(set, itembuffer, states, old_states, sdh, args...)
end