# IronWorkers.jl

Main structure
* db::AbstractDomainBuffer: Specifices a domain to work over, different content if domain is cells, faces, or interfaces
  * Built-in: DomainBuffer, ThreadedDomainBuffer
* w::AbstractWorker
  * Determine what to do at the work_single_x level, where x is determined by `db`
  * Examples: Ferrite.AbstractSparseAssembler, ReAssembler, KeReAssembler, Integrator
  
Typical workflow
At setup, or changes to dh/grid, 
* Create domain buffer, e.g.: `buffer = setup_cell_work(args...)`

At each new timestep/changes to iteration settings
* Update values in `buffer`, e.g.:
   * `set_time_increment!(buffer, Δt)`
   * `update_states!(buffer)` # old = new
For each assembly work or similar
* Create worker, e.g.: `worker = start_assemble(K, r)`
* Do the work: `work!(worker, buffer; a=a, aold=aold)`
Can also have more, e.g. looping over faces, but for std cases, 
this could be handled by the `LoadHandler`.

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
    state # (Updated in element routine)
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

# First dispatch on type of buffer
work_single!(worker, cellbuffer::AbstractCellBuffer) = work_single_cell!(worker, cellbuffer)
work_single!(worker, facebuffer::AbstractFaceBuffer) = work_single_face!(worker, facebuffer)
work_single!(worker, interfacebuffer::AbstractInterfaceBuffer) = work_single_interface!(worker, interfacebuffer)

# Then dispatch on type of worker, i.e.
work_single_cell!(assembler::AbstractAssembler, cb) = assemble_cell!(assembler, cb)
work_single_cell!(integrator::Integrator, cb) = integrate_cell!(integrator, cb)
# etc...

struct DomainBuffer{B,I,S,SDH<:SubDofHandler}
    sdh::SDH
    itembuffer::B  # cell, face, or interface buffer 
    set::Vector{I} # I=Int (cell), I=FaceIndex (face), or
                   # I=NTuple{2,FaceIndex} (interface)
    states::Dict{Int,S}             # To be updated during "work"
    old_states::Dict{Int,S}             # Only for reference
end

struct ThreadedDomainBuffer{B,I,S,SDH<:SubDofHandler}
    sdh::SDH
    itembuffer::TaskLocals{B,B}         # cell, face, or interface buffer 
    chunks::Vector{Vector{Vector{{I}}}} # I=Int (cell), I=FaceIndex (face), or
    set::Vector{I}                      # I=NTuple{2,FaceIndex} (interface)
    states::Dict{Int,S}             # To be updated during "work"
    old_states::Dict{Int,S}             # Only for reference
end
```

Need top level functions like 
`update_states!(::AbstractDomainBuffer)`,
`update_time"(::AbstractDomainBuffer, Δt)`. 
Global vectors still outside? Motivated by them being "global", not belonging
to a certain domain, same as K and r. 


# How to do a Neumann condition with this new setup?
Each new NBC will add one DomainBuffer for each subdofhandler that 
has (a) the field of the given NBC and (2) overlap between the sdh's cellset 
and specified faceset. 

DomainBuffer
* facebuffer: Containing NeumannMaterial and appropriate FaceValues
* set: Standard
* states_*: `Dict{Int,Nothing}()`, i.e. Empty
* sdh # as is


struct NeumannMaterial{F<:Function}
    f::F
    dr::UnitRange # dofrange for given field
end

function face_residual!(re, ae, m::NeumannMaterial, fv, facebuffer)
    for q_point in 1:getnquadpoints(fv)
        dΓ = getdetJdV(fv, q_point)
        x = spatial_coordinate(fv, q_point, getcoordinates(face))
        n = getnormal(fv, q_point)
        b = f(x, time, n)
        for (I, i) in pairs(m.dr)
            δu = shape_value(fv, q_point, i)
            re[i] += (δu ⋅ b) * dΓ
        end
    end
end
