abstract type AbstractDomainBuffer end

get_num_tasks(::AbstractDomainBuffer) = Threads.nthreads() # Default for now

const DomainBuffers = Dict{String,<:AbstractDomainBuffer}
# Acessor functions
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
"""
get_itembuffer(db::DomainBuffers, domain::String) = get_itembuffer(db[domain])

"""
    get_new_state(dbs::Dict{String,AbstractDomainBuffer}, domain::String)
    get_new_state(db::AbstractDomainBuffer)

Get the `states::Dict{Int,S}`, where `S` type of the state for each entity in the domain,
stored in the `db` or `dbs[domain]`.
"""
get_new_state(db::DomainBuffers, domain::String, args...) = get_new_state(db[domain], args...)
"""
    get_old_state(dbs::Dict{String,AbstractDomainBuffer}, domain::String)
    get_old_state(db::AbstractDomainBuffer)

Get the `states::Dict{Int,S}`, where `S` type of the state for each entity in the domain,
stored in the `db` or `dbs[domain]`.
"""
get_old_state(db::DomainBuffers, domain::String, args...) = get_old_state(db[domain], args...)

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

Update the states such that `old_states = new_states` for the states 
stored in `db`
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
    getcellset(dbs::Dict{String,AbstractDomainBuffer}, domain::String)
    getcellset(db::AbstractDomainBuffer)

Get the cellset for `db` or `dbs[domain]`
"""
Ferrite.getcellset(b::DomainBuffers, domain) = getcellset(b[domain])

# Specific AbstractDomainBuffer
# Helpers for constructors
function intersect_set(set::Union{Set{Int},Vector{Int}}, cellset)
    intersected = set===cellset ? set : intersect(set, cellset)
    return sort!(collect(intersected))
end
function intersect_set(set::Union{Set{FaceIndex},Vector{FaceIndex}}, cellset)
    intersected_set = resize!(Vector{FaceIndex}(undef, length(set)), 0)
    for (cellnr, facenr) in set
        cellnr ∈ cellset && push!(intersected_set, FaceIndex(cellnr, facenr))
    end
    return sort!(intersected_set; by=first)
end

struct DomainBuffer{I,B,S,SDH<:SubDofHandler} <: AbstractDomainBuffer
    set::Vector{I}
    itembuffer::B
    new_states::Dict{Int,S} # Always indexed by cell. If desired to have 1 state per e.g. face, need to have e.g. S::Vector{FS}
    old_states::Dict{Int,S} # For interfaces, it is possible/likely that state_here and state_there can be given. 
    sdh::SDH
end
DomainBuffer(set, args...) = DomainBuffer(collect(set), args...)

struct ThreadedDomainBuffer{I,B,S,SDH<:SubDofHandler} <: AbstractDomainBuffer
    chunks::Vector{Vector{Vector{I}}}   # I=Int (cell), I=FaceIndex (face), or
    set::Vector{I}                      # I=NTuple{2,FaceIndex} (interface)
    itembuffer::TaskLocals{B,B}         # cell, face, or interface buffer 
    new_states::Dict{Int,S}             # To be updated during "work"
    old_states::Dict{Int,S}             # Only for reference
    sdh::SDH
end
function ThreadedDomainBuffer(set, itembuffer::AbstractItemBuffer, new_states, old_states, sdh::SubDofHandler, colors_or_chunks=nothing)
    grid = _getgrid(sdh)
    set_vector = collect(set)
    chunks = create_chunks(grid, set_vector, colors_or_chunks)
    itembuffers = TaskLocals(itembuffer)
    return ThreadedDomainBuffer(chunks, set_vector, itembuffers, new_states, old_states, sdh)
end

get_chunks(db::ThreadedDomainBuffer) = db.chunks

const StdDomainBuffer = Union{DomainBuffer, ThreadedDomainBuffer}

Ferrite.getcellset(b::StdDomainBuffer) = b.set
get_dofhandler(b::StdDomainBuffer) = b.sdh.dh
get_itembuffer(b::StdDomainBuffer) = b.itembuffer
get_new_state(b::StdDomainBuffer, cellnum::Int) = fast_getindex(b.new_states, cellnum)
get_old_state(b::StdDomainBuffer, cellnum::Int) = fast_getindex(b.old_states, cellnum)

get_new_state(b::StdDomainBuffer) = b.new_states
get_old_state(b::StdDomainBuffer) = b.old_states
get_material(b::StdDomainBuffer) = get_material(get_base(get_itembuffer(b)))

# Update old_states = new_states after convergence 
function update_states!(b::AbstractDomainBuffer)
    update_states!(get_old_state(b), get_new_state(b))
end

function set_time_increment!(b::StdDomainBuffer, Δt)
    set_time_increment!(get_base(get_itembuffer(b)), Δt)
end

