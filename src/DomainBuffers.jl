abstract type AbstractDomainBuffer end

"""
    DomainBuffers{V<:AbstractDomainBuffer} <: AbstractDict{String,V}

Container for multiple domain buffers. It is publicly fixed (immutable), with 
`getindex`, `haskey`, `get`, `keys`, `values`, `pairs`, `length`, `isempty` and
iteration all forwarded to an internally-wrapped `Dict{String,V}`.

`setindex!` (i.e. `dbs["a"] = db`) and `getproperty` (i.e. `dbs.dict`) are not defined.

Internally, the wrapped `Dict` remains reachable via `_getdict(d)`, and can still be
mutated in place when that is actually safe to do - see `_internal_setindex!`/
`_internal_deletekeys!`.
"""
struct DomainBuffers{V<:AbstractDomainBuffer} <: AbstractDict{String,V}
    dict::Dict{String,V}
    global _unsafe_domainbuffers(d::Dict{String,V}) where {V} = new{V}(d)
end
DomainBuffers(itr) = _unsafe_domainbuffers(Dict(itr))
DomainBuffers{V}() where {V} = _unsafe_domainbuffers(Dict{String,V}())

Base.getproperty(::DomainBuffers, ::Symbol) = error(
    "`DomainBuffers` fields are not accessible via `.`; this type has no public fields")
Base.propertynames(::DomainBuffers) = ()

_getdict(d::DomainBuffers) = getfield(d, :dict)
Base.length(d::DomainBuffers) = length(_getdict(d))
Base.isempty(d::DomainBuffers) = isempty(_getdict(d))
Base.haskey(d::DomainBuffers, k) = haskey(_getdict(d), k)
Base.get(d::DomainBuffers, k, default) = get(_getdict(d), k, default)
Base.getindex(d::DomainBuffers, k) = _getdict(d)[k]
Base.iterate(d::DomainBuffers) = iterate(_getdict(d))
Base.iterate(d::DomainBuffers, state) = iterate(_getdict(d), state)
# `keys`/`values`/`pairs` are intentionally NOT forwarded to the wrapped `Dict`: doing so
# (e.g. `Base.keys(d::DomainBuffers) = keys(_getdict(d))`) would return a `KeySet`/
# `ValueIterator` wrapping the *raw* `Dict`, whose own `.dict` field (unprotected) would hand
# it right back out. Falling through to `AbstractDict`'s generic defaults instead wraps `d`
# itself (so `keys(d).dict`/`values(d).dict` just returns `d`, which is no less protected) -
# at no performance or allocation cost, since those defaults only need `length`/`iterate`.

"""
    _internal_setindex!(d::DomainBuffers, kv::Pair)

Insert `kv` into `d`'s wrapped `Dict` in place (replacing an already-present key - this is a
real `Dict`), bypassing the public immutability.
Usage must be safe wrt. possible multiple references to `d`.
"""
_internal_setindex!(d::DomainBuffers, kv::Pair) = (_getdict(d)[kv.first] = kv.second; d)

"""
    _internal_deletekeys!(d::DomainBuffers, keys)

Remove the entries in `d`'s wrapped `Dict` in place with key in keys,
bypassing the public immutability.
Usage must be safe wrt. possible multiple references to `d`.
"""
_internal_deletekeys!(d::DomainBuffers, keys) = foreach(k -> delete!(_getdict(d), k), keys)

# Accessor functions
"""
    get_dofhandler(dbs::FerriteAssembly.DomainBuffers)
    get_dofhandler(db::AbstractDomainBuffer)
    get_dofhandler(sim::Simulation)

Get the dofhandler stored in `db`. Note that this is the global dofhandler,
and not the `SubDofHandler` that is local to a specific domain.
"""
get_dofhandler(db::DomainBuffers) = get_dofhandler(first(values(db)))

"""
    get_grid(dbs::FerriteAssembly.DomainBuffers)
    get_grid(db::AbstractDomainBuffer)
    get_grid(sim::Simulation)

Get the underlying `grid::AbstractGrid` from the domain buffers
"""
get_grid(db) = Ferrite.get_grid(get_dofhandler(db))

"""
    get_itembuffer(dbs::FerriteAssembly.DomainBuffers, domain::String)
    get_itembuffer(db::AbstractDomainBuffer)
    get_itembuffer(sim::Simulation[, domain::String])

Get the `AbstractItemBuffer` stored in `db` or `dbs[domain]`. 
This internal function might change, but currently the full TaskLocals 
is returned for a ThreadedDomainBuffer (used internally). 
"""
get_itembuffer(db::DomainBuffers, domain::String) = get_itembuffer(db[domain])

"""
    get_state(dbs::FerriteAssembly.DomainBuffers, domain::String)
    get_state(db::Union{AbstractDomainBuffer,FerriteAssembly.DomainBuffers})
    get_state(sim::Simulation[, domain::String])

Get the `states::StateVector{S}`, where `S` is the type of the state for each entity in the
domain, stored in the `db` or `dbs[domain]`. `StateVector{S} <: AbstractDict{Int,S}` is indexed
by cell number and supports the usual read-only `AbstractDict` interface (`keys`, `values`,
`pairs`, iteration, `haskey`, `get`, equality), but is not a `Dict`. If no `domain` is given for
multiple domains, a `Dict{String}` is returned with state variables for each domain.
"""
get_state(db::DomainBuffers, domain::String) = get_state(db[domain])
get_state(db::DomainBuffers) = Dict(key=>get_state(val) for (key,val) in db)

"""
    get_old_state(dbs::FerriteAssembly.DomainBuffers, domain::String)
    get_old_state(db::Union{AbstractDomainBuffer,FerriteAssembly.DomainBuffers})
    get_old_state(sim::Simulation[, domain::String])

Get the `states::StateVector{S}`, where `S` is the type of the state for each entity in the
domain, stored in the `db` or `dbs[domain]`. `StateVector{S} <: AbstractDict{Int,S}` is indexed
by cell number and supports the usual read-only `AbstractDict` interface (`keys`, `values`,
`pairs`, iteration, `haskey`, `get`, equality), but is not a `Dict`. If no `domain` is given for
multiple domains, a `Dict{String}` is returned with state variables for each domain.
"""
get_old_state(db::DomainBuffers, domain::String) = get_old_state(db[domain])
get_old_state(db::DomainBuffers) = Dict(key=>get_old_state(val) for (key,val) in db)

"""
    get_material(dbs::FerriteAssembly.DomainBuffers, domain::String)
    get_material(db::AbstractDomainBuffer)
    get_material(sim::Simulation)

Get the material for the domain represented by `db` or `dbs[domain]`.

**Note:** This always returns the *base* material. If `db` is threaded, `work!` uses
independent copies of this material, so mutating the returned object will not be
reflected when doing threaded assembly. However, in most cases the material is
immutable, and this is not a problem.
"""
get_material(db::DomainBuffers, domain::String) = get_material(db[domain])

# Update functions
"""
    update_states!(db::FerriteAssembly.DomainBuffers; mode::Symbol = :copy)
    update_states!(db::AbstractDomainBuffer; mode::Symbol = :copy)
    update_states!(sim::Simulation; mode::Symbol = :copy)

Update the states such that `old_states = states` (the just-converged values) for the
states stored in `db`.

`mode` selects how this is done:
* `mode = :copy` (default): copies the values from `states` into `old_states`; `states`
  itself is left untouched. This means both `old_states` and `states` correctly hold the
  just-converged values directly after the call — safe to read (e.g. for postprocessing)
  immediately afterwards. If [`create_cell_state`](@ref) returns a *mutable*
  `AbstractArray`, this reuses that array's own storage (no allocation for the array
  itself, though copying non-`isbits` elements into it may still allocate — unless
  [`FerriteAssembly.copy_state!`](@ref) is overloaded for the element type, see below), and
  it must therefore keep the same axes between calls (`ArgumentError` otherwise). Any other
  non-`isbits` cell state must have a [`FerriteAssembly.copy_state`](@ref) or
  [`FerriteAssembly.copy_state!`](@ref) method — otherwise a `MethodError` is thrown. This is
  a **breaking change** from previous releases (which behaved like `mode = :flip`): a
  mutable, non-array cell state without a `copy_state`/`copy_state!` overload that used to
  work now throws; use `mode = :flip` to keep the old behavior for such states.
* `mode = :flip`: cheaply swaps the references of `old_states` and `states` (no copying, no
  allocation, and no `copy_state` requirement — this is the behavior of `update_states!` in
  releases prior to this change). After the call, `states` (the "new" container) holds the
  *stale* values from before this step, not the just-converged ones.
  !!! warning
      Under `mode = :flip`, `states` must not be read again until it has been overwritten by
      the next call to `work!` — including implicitly, e.g. via a default
      [`QuadPointEvaluator`](@ref) reading `get_state`/`s` during postprocessing right after
      `update_states!`. Reading it earlier silently observes the previous step's data. If you
      need to read the just-converged state after updating (e.g. for postprocessing), use the
      default `mode = :copy` instead.
"""
function update_states!(dbs::DomainBuffers; kwargs...)
    for db in values(dbs)
        update_states!(db; kwargs...)
    end
end

"""
    revert_states!(db::FerriteAssembly.DomainBuffers)
    revert_states!(db::AbstractDomainBuffer)
    revert_states!(sim::Simulation)

Update the states such that `states = old_states` for the states stored in `db`,
i.e. the opposite direction of [`update_states!`](@ref). This is useful when 
retrying a time increment after a non-converged solution, when the current (new) 
state is used as an initial guess (typical in staggered solution schemes). Requires
[`FerriteAssembly.copy_state`](@ref) or [`FerriteAssembly.copy_state!`](@ref) for
non `isbits` with the same requirements as stated in
[`update_states!`](@ref) with `mode = :copy`.
"""
function revert_states!(dbs::DomainBuffers)
    for db in values(dbs)
        revert_states!(db)
    end
end

Base.@deprecate set_new_to_old_states! revert_states!

"""
    set_time_increment!(db::FerriteAssembly.DomainBuffers, Δt)
    set_time_increment!(db::AbstractDomainBuffer, Δt)
    set_time_increment!(sim::Simulation, Δt)

Update the time increment stored in `db`, which is passed on to the 
stored `AbstractItemBuffer`
"""
function set_time_increment!(dbs::DomainBuffers, Δt)
    for db in values(dbs)
        set_time_increment!(db, Δt)
    end
end

"""
    replace_material(db::FerriteAssembly.DomainBuffers, replacement_function)
    replace_material(db::AbstractDomainBuffer, replacement_function)

Return a new instance of `db` where as much as possible is copied by reference, and
where the stored material, `m`, is replaced by `replacement_function(m)`.
"""
function replace_material(dbs::DomainBuffers, replacement_function)
    return DomainBuffers(key=>replace_material(db, replacement_function) for (key,db) in dbs)
end

"""
    replace_material(dbs::FerriteAssembly.DomainBuffers, domain::String, replacement_function)

Return a new instance of `dbs` where as much as possible is copied by reference, and
where the material, `m`, of `dbs[domain]` is replaced by `replacement_function(m)`.
Other domains are copied by reference, unchanged.
"""
function replace_material(dbs::DomainBuffers, domain::String, replacement_function)
    haskey(dbs, domain) || throw(ArgumentError("domain \"$domain\" not found in $(collect(keys(dbs)))"))
    return DomainBuffers(key => (key == domain ? replace_material(db, replacement_function) : db) for (key, db) in dbs)
end

"""
    getset(dbs::FerriteAssembly.DomainBuffers, domain::String)
    getset(db::AbstractDomainBuffer)
    getset(sim::Simulation[, domain::String])

Get the set of items stored in `db` or `dbs[domain]`
"""
getset(b::DomainBuffers, domain) = getset(b[domain])

struct DomainBuffer{I,B,SV<:StateVariables,SDH<:SubDofHandler} <: AbstractDomainBuffer
    set::Vector{I}
    itembuffer::B
    states::SV
    sdh::SDH
end

struct ThreadedDomainBuffer{I,B,SV<:StateVariables,SDH<:SubDofHandler} <: AbstractDomainBuffer
    chunks::Vector{Vector{Vector{I}}}   # I=Int (cell), I=FacetIndex (facet), or
    set::Vector{I}                      # I=NTuple{2,FacetIndex} (interface)
    num_tasks::Int
    itembuffer::TaskLocals{B,B}         # cell, facet, or interface buffer
    states::SV
    sdh::SDH
end
function ThreadedDomainBuffer(set, itembuffer::AbstractItemBuffer, states::StateVariables, sdh::SubDofHandler, colors_or_chunks=nothing; num_tasks = Threads.nthreads())
    num_tasks > 0 || throw(ArgumentError("num_tasks must be positive, got $num_tasks"))
    grid = _getgrid(sdh)
    set_vector = collect(set)
    chunks = create_chunks(grid, set_vector, colors_or_chunks)
    itembuffers = TaskLocals(itembuffer; num_tasks)
    return ThreadedDomainBuffer(chunks, set_vector, num_tasks, itembuffers, states, sdh)
end

get_num_tasks(db::ThreadedDomainBuffer) = db.num_tasks
get_num_tasks(dbs::DomainBuffers) = maximum(get_num_tasks, values(dbs))

get_chunks(db::ThreadedDomainBuffer) = db.chunks

const StdDomainBuffer = Union{DomainBuffer, ThreadedDomainBuffer}

#Ferrite.getcellset(b::StdDomainBuffer{Int}) = b.set
#Ferrite.getfacetset(b::StdDomainBuffer{FacetIndex}) = b.set
getset(b::StdDomainBuffer) = b.set

get_dofhandler(b::StdDomainBuffer) = b.sdh.dh
get_itembuffer(b::StdDomainBuffer) = b.itembuffer
get_state(b::StdDomainBuffer, cellnum::Int) = fast_getindex(b.states.new, cellnum)
get_old_state(b::StdDomainBuffer, cellnum::Int) = fast_getindex(b.states.old, cellnum)

get_state(b::StdDomainBuffer) = b.states.new
get_old_state(b::StdDomainBuffer) = b.states.old
get_material(b::StdDomainBuffer) = get_material(get_base(get_itembuffer(b)))

# Update old_states = new_states after convergence
update_states!(b::StdDomainBuffer; kwargs...) = update_states!(b.states; kwargs...)

revert_states!(b::StdDomainBuffer) = revert_states!(b.states)

function set_time_increment!(b::StdDomainBuffer, Δt)
    set_time_increment!(get_base(get_itembuffer(b)), Δt)
end

function replace_material(db::DomainBuffer, replacement_function)
    itembuffer = _replace_material(db.itembuffer, replacement_function)
    return setproperties(db; itembuffer)
end
function replace_material(db::ThreadedDomainBuffer, replacement_function)
    base_ibuf = _replace_material(get_base(db.itembuffer), replacement_function)
    task_ibuf = map(ibuf->_replace_material(ibuf, replacement_function), get_locals(db.itembuffer))
    return setproperties(db; itembuffer = TaskLocals(base_ibuf, task_ibuf))
end

# Experimental: Insert new states, allows reusing the buffer for multiple simulations with same
# initial state (grid, dh, etc.), but which experience different loading. Typically for RVE simulations. 
function replace_states!(dbs::DomainBuffers, states::Dict{String, <:StateVariables})
    keys(dbs) == keys(states) || throw(ArgumentError("keys of dictionaries don't match"))
    for (key, db) in dbs
        replace_states!(db, states[key])
    end
end

replace_states!(db::StdDomainBuffer, states::StateVariables) = replace_states!(db.states, states)
