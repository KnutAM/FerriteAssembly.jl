# IronWorkers.jl
Naming
* item -> entity?
  
  
## ItemBuffer
```julia
struct CellBuffer <: AbstractItemBuffer
    # (Normally) fixed during simulation (set at start)
    dofrange
    material
    # Set at beginning of loop
    Δt
    user_data
    # Re-point before element routine
    new_state # (Updated in element routine)
    old_state 
    # Reinit before element routine
    ae_old
    ae
    dofs
    coords
    cellvalues
    cellid
    # Reset before element routine (zero)
    Ke
    re
    # Do nothing (may be changed in element routine)
    user_cache
end

update_time!(c::CellBuffer, Δt)
update_material!(c::CellBuffer, m)
```
## General nested levels 

```julia
function work!(worker, buffers::Dict{String,DomainBuffer}; kwargs...)
    for buffer in values(buffers)
        work_domain_sequential!(worker, buffer; kwargs...)
    end
end

function work_domain_sequential!(worker, domainbuffer; kwargs...)
    itembuffer = get_itembuffer(domainbuffer)
    for cellnr in getcellset(domainbuffer)
        reinit!(itembuffer, cellnr, domainbuffer; kwargs...)
        work_single!(worker, itembuffer)
    end
end

function work!(worker, buffers::Dict{String,ThreadedDomainBuffer}; kwargs...)
    workers = TaskLocals(worker)
    for buffer in values(buffers)
        work_domain_threaded!(workers, buffer; kwargs...)
    end
end
function work_domain_threaded!(workers, domainbuffer; kwargs...)
    itembuffers = get_itembuffer(domainbuffer) #::TaskLocals
    scatter!(itembuffers)
    scatter!(workers)
    num_tasks = get_num_tasks(domainbuffer) # Default to Threads.nthreads()
    for chunk_vector in get_chunk_vectors(domainbuffer)
        taskchunks = TaskChunks(chunk_vector)
        Base.Experimental.@sync begin 
            for taskid in 1:num_tasks
                itembuffer = get_local(itembuffers, taskid)
                worker = get_local(workers, taskid)
                Threads.@spawn begin
                    while true
                        taskchunk = get_chunk(taskchunks)
                        isempty(taskchunk) && break
                        for cellnr in taskchunk
                            reinit!(itembuffer, cellnr, domainbuffer; kwargs...)
                            work_single!(worker, itembuffer)
                        end
                    end
                end #spawn
            end #taskid
        end #sync
    end #colors
    gather!(cellbuffers)
    gather!(assemblers)
end

# Dispatch on type of buffer
work_single!(worker, cellbuffer::AbstractCellBuffer) = work_single_cell!(worker, cellbuffer)
work_single!(worker, facebuffer::AbstractFaceBuffer) = work_single_face!(worker, facebuffer)
work_single!(worker, interfacebuffer::AbstractInterfaceBuffer) = work_single_interface!(worker, interfacebuffer)

# Then dispatch on type of worker, i.e.
work_single_cell!(assembler::AbstractAssembler, cb) = ...

struct DomainBuffer{B,I,S,SDH<:SubDofHandler}
    sdh::SDH
    itembuffer::B  # cell, face, or interface buffer 
    set::Vector{I} # I=Int (cell), I=FaceIndex (face), or
                   # I=NTuple{2,FaceIndex} (interface)
    new_states::Dict{Int,S}             # To be updated during "work"
    old_states::Dict{Int,S}             # Only for reference
end

struct ThreadedDomainBuffer{B,I,S,SDH<:SubDofHandler}
    sdh::SDH
    itembuffer::TaskLocals{B,B}         # cell, face, or interface buffer 
    chunks::Vector{Vector{Vector{{I}}}} # I=Int (cell), I=FaceIndex (face), or
    set::Vector{I}                      # I=NTuple{2,FaceIndex} (interface)
    new_states::Dict{Int,S}             # To be updated during "work"
    old_states::Dict{Int,S}             # Only for reference
end
```

Need top level functions like 
`update_states!(::AbstractDomainBuffer)`,
`update_time"(::AbstractDomainBuffer, Δt)`. 
Global vectors still outside? Motivated by them being "global", not belonging
to a certain domain, same as K and r. 