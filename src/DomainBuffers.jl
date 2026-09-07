abstract type AbstractDomainBuffer end

const DomainBuffers = Dict{String, <:AbstractDomainBuffer}

# Accessor functions
"""
    get_dofhandler(dbs::Dict{String,AbstractDomainBuffer})
    get_dofhandler(db::AbstractDomainBuffer)
    get_dofhandler(sim::Simulation)

Get the dofhandler stored in `db`. Note that this is the global dofhandler,
and not the `SubDofHandler` that is local to a specific domain.
"""
get_dofhandler(db::DomainBuffers) = get_dofhandler(first(values(db)))

"""
    get_grid(dbs::Dict{String,AbstractDomainBuffer})
    get_grid(db::AbstractDomainBuffer)
    get_grid(sim::Simulation)

Get the underlying `grid::AbstractGrid` from the domain buffers
"""
get_grid(db) = Ferrite.get_grid(get_dofhandler(db))

"""
    get_itembuffer(dbs::Dict{String,AbstractDomainBuffer}, domain::String)
    get_itembuffer(db::AbstractDomainBuffer)
    get_itembuffer(sim::Simulation[, domain::String])

Get the `AbstractItemBuffer` stored in `db` or `dbs[domain]`. 
This internal function might change, but currently the full TaskLocals 
is returned for a ThreadedDomainBuffer (used internally). 
"""
get_itembuffer(db::DomainBuffers, domain::String) = get_itembuffer(db[domain])

"""
    get_state(dbs::Dict{String,AbstractDomainBuffer}, domain::String)
    get_state(db::Union{AbstractDomainBuffer,Dict{String,AbstractDomainBuffer}})
    get_state(sim::Simulation[, domain::String])

Get the `states::Dict{Int,S}`, where `S` type of the state for each entity in the domain,
stored in the `db` or `dbs[domain]`. If no `domain` is given for multiple domains, a 
Dict{String} is returned with state variables for each domain
"""
get_state(db::DomainBuffers, domain::String) = get_state(db[domain])
get_state(db::DomainBuffers) = Dict(key=>get_state(val) for (key,val) in db)

"""
    get_old_state(dbs::Dict{String,AbstractDomainBuffer}, domain::String)
    get_old_state(db::Union{AbstractDomainBuffer,Dict{String,AbstractDomainBuffer}})
    get_old_state(sim::Simulation[, domain::String])

Get the `states::Dict{Int,S}`, where `S` type of the state for each entity in the domain,
stored in the `db` or `dbs[domain]`. If no `domain` is given for multiple domains, a 
Dict{String} is returned with state variables for each domain
"""
get_old_state(db::DomainBuffers, domain::String) = get_old_state(db[domain])
get_old_state(db::DomainBuffers) = Dict(key=>get_old_state(val) for (key,val) in db)

"""
    get_material(dbs::Dict{String,AbstractDomainBuffer}, domain::String)
    get_material(db::AbstractDomainBuffer)
    get_material(sim::Simulation)

Get the material for the domain represented by `db` or `dbs[domain]`.
"""
get_material(db::DomainBuffers, domain::String) = get_material(db[domain])

# Update functions
"""
    update_states!(db::Dict{String,AbstractDomainBuffer}; mode::Symbol = :copy)
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
  itself, though copying non-`isbits` elements into it may still allocate via
  `copy_state`), and it must therefore keep the same axes between calls (`ArgumentError`
  otherwise). Any other non-`isbits` cell state must have a
  [`FerriteAssembly.copy_state`](@ref) method — otherwise a `MethodError` is thrown. This is
  a **breaking change** from previous releases (which behaved like `mode = :flip`): a
  mutable, non-array cell state without a `copy_state` overload that used to work now
  throws; use `mode = :flip` to keep the old behavior for such states.
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
    set_new_to_old_states!(db::Dict{String,AbstractDomainBuffer})
    set_new_to_old_states!(db::AbstractDomainBuffer)
    set_new_to_old_states!(sim::Simulation)

!!! warning "Deprecated"
    `set_new_to_old_states!` is deprecated and may be removed in a future release.

Update the states such that `states = old_states` for the states stored in `db`,
i.e. the opposite direction of [`update_states!`](@ref). This is useful for
resetting the current (new) state to the last converged (old) state, e.g. when
retrying a time increment after a non-converged solution, without having to
reassemble.

Unlike `update_states!(db; mode=:flip)`, this method does not swap references between
`old_states` and `states`, but copies values from `old_states` into the existing `states`
containers. If [`create_cell_state`](@ref) returns a *mutable* `AbstractArray`
(`ismutable(state) == true`), each element is copied individually; otherwise (including an
immutable `AbstractArray`) the whole cell state is copied. In both cases, values for which
`isbits(value) == true` are copied by identity (no allocation); any other value is copied
with [`FerriteAssembly.copy_state`](@ref), which has no default method and must be
overloaded for that value's type — otherwise a `MethodError` is thrown.

!!! note
    When [`create_cell_state`](@ref) returns a mutable `AbstractArray`, that array in
    `states` is updated in-place and must therefore have the same axes as the corresponding
    array in `old_states` — an `ArgumentError` is thrown otherwise. An immutable
    `AbstractArray` cell state (e.g. built from `NTuple`s or `StaticArrays`) is instead
    replaced wholesale, like any other non-array
    cell state.
"""
function set_new_to_old_states!(dbs::DomainBuffers)
    for db in values(dbs)
        set_new_to_old_states!(db)
    end
end

"""
    set_time_increment!(db::Dict{String,AbstractDomainBuffer}, Δt)
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
    replace_material(db::Dict{String,AbstractDomainBuffer}, replacement_function)
    replace_material(db::AbstractDomainBuffer, replacement_function)

Return a new instance of `db` where as much as possible is copied by reference, and 
where the stored material, `m`, is replaced by `replacement_function(m)`.
"""
function replace_material(dbs::DomainBuffers, replacement_function)
    return Dict(key=>replace_material(db, replacement_function) for (key,db) in dbs)
end

"""
    couple_buffers(dbs::Dict{String, <:AbstractDomainBuffer}; kwargs::Dict{String, <:AbstractDomainBuffer}...)
    couple_buffers(db::AbstractDomainBuffer; kwargs::AbstractDomainBuffer...)

Return new buffer(s) that are coupled with the buffers provided as keyword arguments. The key is used in 
[`get_coupled_buffer`](@ref) to get the coupled itembuffer, such that its values may be queried. 

!!! note
    This functionality assumes that each setup has the same grid, and in case of multiple domains, these should also 
    match.
"""
function couple_buffers(dbs::DomainBuffers; kwargs...)
    return Dict(
        key => (all(haskey(v, key) for (_, v) in kwargs) ? 
            couple_buffers(db; (k => v[key] for (k, v) in kwargs)...) :
            db) for (key, db) in dbs)
    #return Dict(key => couple_buffers(db; (k => v[key] for (k, v) in kwargs)...) for (key, db) in dbs)
end

"""
    getset(dbs::Dict{String,AbstractDomainBuffer}, domain::String)
    getset(db::AbstractDomainBuffer)
    getset(sim::Simulation[, domain::String])

Get the set of items stored in `db` or `dbs[domain]`
"""
getset(b::DomainBuffers, domain) = getset(b[domain])

struct DomainBuffer{I,B,S,SDH<:SubDofHandler} <: AbstractDomainBuffer
    set::Vector{I}
    itembuffer::B
    states::StateVariables{S}
    sdh::SDH
end

struct ThreadedDomainBuffer{I,B,S,SDH<:SubDofHandler} <: AbstractDomainBuffer
    chunks::Vector{Vector{Vector{I}}}   # I=Int (cell), I=FacetIndex (facet), or
    set::Vector{I}                      # I=NTuple{2,FacetIndex} (interface)
    num_tasks::Int
    itembuffer::TaskLocals{B,B}         # cell, facet, or interface buffer 
    states::StateVariables{S}
    sdh::SDH
end
function ThreadedDomainBuffer(set, itembuffer::AbstractItemBuffer, states::StateVariables, sdh::SubDofHandler, colors_or_chunks=nothing; num_tasks = Threads.nthreads())
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

set_new_to_old_states!(b::StdDomainBuffer) = set_new_to_old_states!(b.states)

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

function couple_buffers(db::DomainBuffer; kwargs...)
    itembuffer = couple_buffers(db.itembuffer; (k => v.itembuffer for (k, v) in kwargs)...)
    return setproperties(db; itembuffer)
end

function couple_buffers(db::ThreadedDomainBuffer; kwargs...)
    base_ibuf = couple_buffers(get_base(db.itembuffer); (k => get_base(v.itembuffer) for (k, v) in kwargs)...)
    task_ibuf = map(enumerate(get_locals(db.itembuffer))) do (i, ibuf)
        couple_buffers(ibuf; (k => get_local(v.itembuffer, i) for (k, v) in kwargs)...)
    end
    return setproperties(db; itembuffer = TaskLocals(base_ibuf, task_ibuf))
end

# Experimental: Insert new states, allows reusing the buffer for multiple simulations with same 
# initial state (grid, dh, etc.), but which experience different loading. Typically for RVE simulations. 
function replace_states!(dbs::Dict{String, <:AbstractDomainBuffer}, states::Dict{String, <:StateVariables})
    keys(dbs) == keys(states) || throw(ArgumentError("keys of dictionaries don't match"))
    for (key, db) in dbs
        replace_states!(db, states[key])
    end
end

replace_states!(db::StdDomainBuffer, states::StateVariables) = replace_states!(db.states, states)
