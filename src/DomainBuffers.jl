abstract type AbstractDomainBuffer end
# Acessor functions
function get_dofhandler end
function get_itembuffer end 
function get_new_state end 
function get_old_state end 

# Update functions 
function update_states! end
function update_time_increment! end 

struct DomainBuffer{B,I,S,SDH<:SubDofHandler} <: AbstractDomainBuffer
    itembuffer::B
    set::Vector{I}
    new_states::Dict{Int,S} # Always indexed by cell. If desired to have 1 state per e.g. face, need to have e.g. S::Vector{FS}
    old_states::Dict{Int,S} # For interfaces, it is possible/likely that state_here and state_there can be given. 
    sdh::SDH
end

struct ThreadedDomainBuffer{B,I,S,SDH<:SubDofHandler} <: AbstractDomainBuffer
    itembuffer::TaskLocals{B,B}         # cell, face, or interface buffer 
    chunks::Vector{Vector{Vector{{I}}}} # I=Int (cell), I=FaceIndex (face), or
    set::Vector{I}                      # I=NTuple{2,FaceIndex} (interface)
    new_states::Dict{Int,S}             # To be updated during "work"
    old_states::Dict{Int,S}             # Only for reference
    sdh::SDH
end

const StdDomainBuffer = Union{DomainBuffer, ThreadedDomainBuffer}

get_dofhandler(b::StdDomainBuffer) = b.sdh.dh
get_itembuffer(b::StdDomainBuffer) = b.itembuffer
get_new_state(b::StdDomainBuffer, cellnum::Int) = fast_getindex(b.new_states, cellnum)
get_old_state(b::StdDomainBuffer, cellnum::Int) = fast_getindex(b.old_states, cellnum)

get_new_state(b::StdDomainBuffer) = b.new_states
get_old_state(b::StdDomainBuffer) = b.old_states

# Update old_states = new_states after convergence 
function update_states!(b::AbstractDomainBuffer)
    update_states!(get_old_state(b), get_new_state(b))
end
function update_states!(bs::Dict{String,<:AbstractDomainBuffer})
    foreach(update_states!, values(bs))
end

function update_time_increment!(bs::Dict{String,<:AbstractDomainBuffer}, Δt)
    for b in values(bs)
        update_time_increment!(b, Δt)
    end
end
function update_time_increment!(b::DomainBuffer, Δt)
    update_time_increment!(get_itembuffer(b), Δt)
end
function update_time_increment!(b::ThreadedDomainBuffer, Δt)
    update_time_increment!(get_base(get_itembuffer(b)), Δt)
end

