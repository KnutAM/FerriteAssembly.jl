"""
    DomainSpec(sdh::SubDofHandler, material, fe_values; set=getcellset(sdh), colors_or_chunks=nothing, user_data=nothing)
    DomainSpec(dh::DofHandler, material, fe_values; set=1:getncells(dh), colors=nothing, chunks=nothing, user_data=nothing)

Create a `DomainSpec` that can be used to set up a domain buffer. 
The element type of `set` determines the type of domain 
* `Int`:       cell domain, `fe_values` should be `CellValues` (or a `NamedTuple` with `CellValues` as elements)
* `FaceIndex`: face domain, `fe_values` should be `FaceValues` (or a `NamedTuple` with `FaceValues` as elements)

`material` is used for dispatch on the utilized `worker`'s function. 
The `user_data` is just passed along (by reference) everywhere, and is accessible 
via the `ItemBuffer` (e.g. `CellBuffer`) given to the `worker`'s function.

For multithreading, `colors::Union{Vector{Vector{Int}},Vector{Vector{eltype(set)}}` **or** `chunks::Vector{Vector{Vector{eltype(set)}}}` 
may be passed in order to partition the grid. Internally, if not given, `chunks` will be created from a colored grid (either given here,
or calculated using Ferrite's default coloring algorithm). Note that the chunks must refer to the values in the set (exactly same), whereas
colors can be for the entire grid, although it is usually better to color each cellset individually.
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
    return DomainSpec(SubDofHandler(dh), args...; kwargs...)
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