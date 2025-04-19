abstract type AbstractDomainBuffer end

get_num_tasks(::AbstractDomainBuffer) = Threads.nthreads() # Default for now

const DomainBuffers = Dict{String,<:AbstractDomainBuffer}

# Accessor functions
"""
    get_dofhandler(dbs::Dict{String,AbstractDomainBuffer})
    get_dofhandler(db::AbstractDomainBuffer)

Get the dofhandler stored in `db`. Note that this is the global dofhandler,
and not the `SubDofHandler` that is local to a specific domain.
"""
get_dofhandler(db::DomainBuffers) = get_dofhandler(first(values(db)))

"""
    get_itembuffer(dbs::Dict{String,AbstractDomainBuffer}, domain::String)
    get_itembuffer(db::AbstractDomainBuffer)

Get the `AbstractItemBuffer` stored in `db` or `dbs[domain]`. 
This internal function might change, but currently the full TaskLocals 
is returned for a ThreadedDomainBuffer (used internally). 
"""
get_itembuffer(db::DomainBuffers, domain::String) = get_itembuffer(db[domain])

"""
    get_state(dbs::Dict{String,AbstractDomainBuffer}, domain::String)
    get_state(db::Union{AbstractDomainBuffer,Dict{String,AbstractDomainBuffer}})

Get the `states::Dict{Int,S}`, where `S` type of the state for each entity in the domain,
stored in the `db` or `dbs[domain]`. If no `domain` is given for multiple domains, a 
Dict{String} is returned with state variables for each domain
"""
get_state(db::DomainBuffers, domain::String) = get_state(db[domain])
get_state(db::DomainBuffers) = Dict(key=>get_state(val) for (key,val) in db)

"""
    get_old_state(dbs::Dict{String,AbstractDomainBuffer}, domain::String)
    get_old_state(db::Union{AbstractDomainBuffer,Dict{String,AbstractDomainBuffer}})

Get the `states::Dict{Int,S}`, where `S` type of the state for each entity in the domain,
stored in the `db` or `dbs[domain]`. If no `domain` is given for multiple domains, a 
Dict{String} is returned with state variables for each domain
"""
get_old_state(db::DomainBuffers, domain::String) = get_old_state(db[domain])
get_old_state(db::DomainBuffers) = Dict(key=>get_old_state(val) for (key,val) in db)

"""
    get_material(dbs::Dict{String,AbstractDomainBuffer}, domain::String)
    get_material(db::AbstractDomainBuffer)

Get the material for the domain represented by `db` or `dbs[domain]`.
"""
get_material(db::DomainBuffers, domain::String) = get_material(db[domain])

# Update functions
"""
    update_states!(db::Dict{String,AbstractDomainBuffer})
    update_states!(db::AbstractDomainBuffer)

Update the states such that `old_states = states` for the states 
stored in `db`.

This method tries to avoid allocating new values where possible. 
Currently, if [`create_cell_state`](@ref) returns `T` or `Vector{T}` where `isbitstype(T)`, this works.
If needed/wanted, it should be relatively easy to provide an interface to make it possible to have allocation free 
for custom cell states.
"""
function update_states!(dbs::DomainBuffers)
    for db in values(dbs)
        update_states!(db)
    end
end

"""
    set_time_increment!(db::Dict{String,AbstractDomainBuffer}, Δt)
    set_time_increment!(db::AbstractDomainBuffer, Δt)

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

Couple the buffer(s) with one or more other buffers provided as keyword arguments. The key can then be used 
to get the itembuffer inside the element routines from the current itembuffer, allowing it to be updated.

!!! note 
    This functionality assumes that each setup has the same grid, and in case of multiple domains, these should also 
    match.
"""
function couple_buffers(dbs::DomainBuffers; kwargs...)
    return Dict(key => couple_buffers(db; (k => v[key] for (k, v) in kwargs)...) for (key, db) in dbs)
end

"""
    getset(dbs::Dict{String,AbstractDomainBuffer}, domain::String)
    getset(db::AbstractDomainBuffer)

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
    itembuffer::TaskLocals{B,B}         # cell, facet, or interface buffer 
    states::StateVariables{S}
    sdh::SDH
end
function ThreadedDomainBuffer(set, itembuffer::AbstractItemBuffer, states::StateVariables, sdh::SubDofHandler, colors_or_chunks=nothing)
    grid = _getgrid(sdh)
    set_vector = collect(set)
    chunks = create_chunks(grid, set_vector, colors_or_chunks)
    itembuffers = TaskLocals(itembuffer)
    return ThreadedDomainBuffer(chunks, set_vector, itembuffers, states, sdh)
end

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
update_states!(b::StdDomainBuffer) = update_states!(b.states)

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
