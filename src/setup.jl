"""
    GridDomain(sdh::SubDofHandler, material, fe_values; set=getcellset(sdh), colors_or_chunks=nothing, user_data=nothing)
    GridDomain(dh::DofHandler, material, fe_values; set=1:getncells(dh), colors=nothing, chunks=nothing, user_data=nothing)

Create a `GridDomain` that can be used to set up a domain buffer. 
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
struct GridDomain{I}
    sdh::SubDofHandler
    material::Any
    fe_values::Any
    set::Union{AbstractVector{I},AbstractSet{I}}
    colors_or_chunks::Any
    user_data::Any
end
function GridDomain(sdh::SubDofHandler, material, values; set=getcellset(sdh), colors=nothing, chunks=nothing, user_data=nothing)
    intersected_set = set === getcellset(sdh) ? set : intersect_with_cellset(set, getcellset(sdh))
    colors_or_chunks = chunks!==nothing ? chunks : colors
    return GridDomain(sdh, material, values, intersected_set, colors_or_chunks, user_data)
end
function GridDomain(dh::DofHandler, args...; kwargs...)
    return GridDomain(SubDofHandler(dh), args...; kwargs...)
end

intersect_with_cellset(set::Union{AbstractSet{Int},AbstractVector{Int}}, cellset) = intersect(set, cellset)


"""
    setup_domainbuffers(gds::Dict{String,GridDomain}; kwargs...)

Setup multiple domain buffers, one for each `GridDomain` in `gds`.
See [`setup_domainbuffer`](@ref) for description of the keyword arguments.
"""
function setup_domainbuffers(gds::Dict{String,<:GridDomain}; kwargs...)
    return Dict(name => setup_domainbuffer(gd; kwargs...) for (name, gd) in gds)
end

"""
    setup_domainbuffer(gd::GridDomain; anew=nothing, threading=false, autodiffbuffer=false)

Setup a domain buffer for a single grid domain, `gd`. 
* `anew::Vector`: The global degree of freedom values are used to pass the 
  local element dof values to the [`create_cell_state`](@ref) function, 
  making it possible to create the initial state dependent on the initial 
  conditions for the field variables. 
* `threading`: Should a `ThreadedDomainBuffer` be created to work the grid multithreaded
  if supported by the used `worker`?
* `autodiffbuffer`: Should a custom itembuffer be used to speed up the automatic 
  differentiation (if supported by the itembuffer)
"""
function setup_domainbuffer(gd::GridDomain; threading=Val(false), kwargs...)
    thrd = isa(threading, Val) ? threading : Val(threading)
    return _setup_domainbuffer(thrd, gd; kwargs...)
end

create_states(gd::GridDomain{Int}, a) = create_states(gd.sdh, gd.material, gd.fe_values, a, gd.set, create_dofrange(gd.sdh))
create_states(::GridDomain{FaceIndex}, ::Any) = Dict{Int,Nothing}()

function setup_itembuffer(adb, gd::GridDomain{FaceIndex}, args...)
    dofrange = create_dofrange(gd.sdh)
    return setup_facebuffer(adb, gd.sdh, gd.fe_values, gd.material, dofrange, gd.user_data)
end
function setup_itembuffer(adb, gd::GridDomain{Int}, states)
    dofrange = create_dofrange(gd.sdh)
    return setup_cellbuffer(adb, gd.sdh, gd.fe_values, gd.material, first(values(states)), dofrange, gd.user_data)
end

function _setup_domainbuffer(threaded, gd; anew=nothing, autodiffbuffer=Val(false))
    new_states = create_states(gd, anew)
    old_states = create_states(gd, anew)
    itembuffer = setup_itembuffer(autodiffbuffer, gd, new_states)
    return _setup_domainbuffer(threaded, gd.set, itembuffer, new_states, old_states, gd.sdh, gd.colors_or_chunks)
end

# Type-unstable switch
function _setup_domainbuffer(threaded::Bool, args...)
    return _setup_domainbuffer(Val(threaded), args...)
end
# Sequential
function _setup_domainbuffer(::Val{false}, set, itembuffer, new_states, old_states, sdh, args...)
    return DomainBuffer(set, itembuffer, new_states, old_states, sdh)
end
# Threaded
function _setup_domainbuffer(::Val{true}, set, itembuffer, new_states, old_states, sdh, args...)
    return ThreadedDomainBuffer(set, itembuffer, new_states, old_states, sdh, args...)
end