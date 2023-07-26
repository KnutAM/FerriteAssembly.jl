
function work!(worker, buffers::Dict{String,DomainBuffer}; kwargs...)
    for buffer in values(buffers)
        work_domain_sequential!(worker, buffer; kwargs...)
    end
end
function work!(worker, buffer::DomainBuffer; kwargs...)
    work_domain_sequential!(worker, buffer; kwargs...)
end
function work!(worker, buffers::Dict{String,ThreadedDomainBuffer}; kwargs...)
    if can_thread(worker)
        workers = TaskLocals(worker)
        for buffer in values(buffers)
            work_domain_threaded!(workers, buffer; kwargs...)
        end
    else
        for buffer in values(buffers)
            work_domain_sequential!(worker, buffer; kwargs...)
        end
    end
end
function work!(worker, buffer::ThreadedDomainBuffer; kwargs...)
    if can_thread(worker)
        workers = TaskLocals(worker)
        work_domain_threaded!(workers, buffer; kwargs...)
    else
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
                        taskchunk = get_chunk(taskchunks) # Vector{Int}
                        isempty(taskchunk) && break
                        for cellnr in taskchunk
                            reinit!(itembuffer, domainbuffer, cellnr; kwargs...)
                            work_single!(worker, itembuffer)
                        end # cellnr
                    end #chunk
                end #spawn
            end #taskid
        end #sync
    end #chunk_vectors
    gather!(itembuffers)
    gather!(workers)
end


"""
    doassemble!(
        assembler, new_states::Dict{Int}, buffer::AbstractDomainBuffer; 
        a=nothing, aold=nothing, old_states=nothing, Δt=NaN)

Use `assembler` to assemble a single domain described by `buffer`, and update `new_states` if dictated by the assembler. 

    doassemble!(
        assembler, new_states::Dict{String}, buffers::Dict{String,AbstractDomainBuffer}; 
        a=nothing, aold=nothing, old_states=nothing, Δt=NaN)

Use `assembler` to assemble the domains described by `buffers`, and update `new_states` if dictated by the assembler.

The keyword arguments work as follows:
- `a, aold`: Global degree of freedom vectors. If `nothing`, `NaN` values are passed to the element routine.
- `old_states`: Old state variables. If `nothing`, [`get_old_state`](@ref FerriteAssembly.get_old_state) 
  will always return the initial state. 
- `Δt`: The value returned by [`get_time_increment`](@ref FerriteAssembly.get_time_increment)
"""
function doassemble!(assembler, new_states, buffer::DomainBuffer; kwargs...)
    sequential_assemble_domain!(assembler, new_states, buffer; kwargs...)
end
function doassemble!(assembler, new_states, buffer::ThreadedDomainBuffer; kwargs...)
    if can_thread(assembler)
        threaded_assemble_domain!(TaskLocals(assembler), new_states, buffer; kwargs...)
    else
        sequential_assemble_domain!(assembler, new_states, buffer; kwargs...)
    end
end

# Assemble multiple domains
function doassemble!(assembler, new_states, buffers::Dict{String,<:DomainBuffer}; kwargs...)
    sequential_assemble_domains!(assembler, new_states, buffers; kwargs...)
end

function doassemble!(assembler, new_states, buffers::Dict{String,<:ThreadedDomainBuffer}; kwargs...)
    if can_thread(assembler)
        threaded_assemble_domains!(TaskLocals(assembler), new_states, buffers; kwargs...)
    else
        sequential_assemble_domains!(assembler, new_states, buffers; kwargs...)
    end
end
function sequential_assemble_domains!(assembler, new_states, buffers::Dict{String}; old_states=nothing, kwargs...)
    for key in keys(new_states)
        skip_this_domain(assembler, key) && continue
        domain_old_states = isnothing(old_states) ? nothing : fast_getindex(old_states, key)
        sequential_assemble_domain!(assembler, new_states[key], buffers[key]; old_states=domain_old_states, kwargs...)
    end
end
function threaded_assemble_domains!(assembler::TaskLocals, new_states, buffers::Dict{String,<:ThreadedDomainBuffer}; old_states=nothing, kwargs...)
    for key in keys(new_states)
        skip_this_domain(get_base(assembler), key) && continue
        domain_old_states = isnothing(old_states) ? nothing : fast_getindex(old_states, key)
        threaded_assemble_domain!(assembler, new_states[key], buffers[key]; old_states=domain_old_states, kwargs...)
    end
end

function sequential_assemble_domain!(assembler, new_states, buffer; a=nothing, aold=nothing, old_states=nothing, Δt=NaN)
    cellbuffer = isa(buffer, DomainBuffer) ? buffer.cellbuffer : get_base(buffer.cellbuffers)
    update_time!(cellbuffer, Δt)
    for cellnr in buffer.cellset
        assemble_cell!(assembler, new_states, cellbuffer, buffer.sdh, cellnr, a, aold, old_states)
    end
end

function threaded_assemble_domain!(assemblers::TaskLocals, new_states, buffer::ThreadedDomainBuffer; a=nothing, aold=nothing, old_states=nothing, Δt=NaN)
    cellbuffers = buffer.cellbuffers #::TaskLocals
    sdh = buffer.sdh # SubDofHandler
    update_time!(get_base(cellbuffers), Δt)
    scatter!(cellbuffers)
    scatter!(assemblers)
    num_tasks = Threads.nthreads()
    for saved_set_chunks in buffer.saved_set_chunks
        set_chunks = TaskChunks(saved_set_chunks)
        Base.Experimental.@sync begin 
            for taskid in 1:num_tasks
                cellbuffer = get_local(cellbuffers, taskid)
                assembler = get_local(assemblers, taskid)
                Threads.@spawn begin
                    while true
                        set_chunk = get_chunk(set_chunks)
                        isempty(set_chunk) && break
                        for cellnr in set_chunk
                            assemble_cell!(assembler, new_states, cellbuffer, sdh, cellnr, a, aold, old_states)
                        end
                    end
                end #spawn
            end #taskid
        end #sync
    end #colors
    gather!(cellbuffers)
    gather!(assemblers)
end

"""
    assemble_cell!(assembler, new_states, cellbuffer, sdh::SubDofHandler, cellnr, a, aold, old_states)

Internal function to that reinitializes the `cellbuffer` and calls [`assemble_cell_reinited!`](@ref).
"""
function assemble_cell!(assembler, new_states, cellbuffer, sdh::SubDofHandler, cellnr, a, aold, old_states)
    reinit!(cellbuffer, sdh.dh, cellnr, a, aold, old_states)
    cell_state = fast_getindex(new_states, cellnr)
    assemble_cell_reinited!(assembler, cell_state, cellbuffer)
end