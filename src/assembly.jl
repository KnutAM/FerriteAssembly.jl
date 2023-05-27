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

threaded_assemble_domain!(args...; kwargs...) = threaded_assemble_domain4!(args...; kwargs...)

# Current (does not support threaded due to :static scheduling)
function threaded_assemble_domain1!(assembler::TaskLocals, new_states, buffer::ThreadedDomainBuffer; a=nothing, aold=nothing, old_states=nothing, Δt=NaN)
    cellbuffers = buffer.cellbuffers #::TaskLocals
    update_time!(get_base(cellbuffers), Δt)
    scatter!(cellbuffers)
    scatter!(assembler)
    for cellset in buffer.colors
        Threads.@threads :static for cellnr in cellset
            id = Threads.threadid()
            assemble_cell!(get_local(assembler, id), new_states, get_local(cellbuffers, id), buffer.sdh, cellnr, a, aold, old_states)
        end
    end
    gather!(cellbuffers)
    gather!(assembler)
end

# Channel implementation
function threaded_assemble_domain2!(assemblers::TaskLocals, new_states, buffer::ThreadedDomainBuffer; a=nothing, aold=nothing, old_states=nothing, Δt=NaN)
    cellbuffers = buffer.cellbuffers #::TaskLocals
    sdh = buffer.sdh # SubDofHandler
    update_time!(get_base(cellbuffers), Δt)
    scatter!(cellbuffers)
    scatter!(assemblers)
    num_tasks = Threads.nthreads()
    queue_length = typemax(Int)# 2*num_tasks
    for saved_set_chunks in buffer.saved_set_chunks
        set_chunks = Channel{Vector{Int}}(queue_length)
        # Base.Experimental.@sync exits immediately on exception (allow Ctrl+C)
        Base.Experimental.@sync begin 
            Threads.@spawn begin
                for set_chunk in saved_set_chunks
                    put!(set_chunks, set_chunk)
                end
                close(set_chunks)
            end 
            for taskid in 1:num_tasks
                cellbuffer = get_local(cellbuffers, taskid)
                assembler = get_local(assemblers, taskid)
                Threads.@spawn begin
                    for set_chunk in set_chunks
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

function threaded_assemble_domain3!(assemblers::TaskLocals, new_states, buffer::ThreadedDomainBuffer; a=nothing, aold=nothing, old_states=nothing, Δt=NaN)
    cellbuffers = buffer.cellbuffers #::TaskLocals
    sdh = buffer.sdh # SubDofHandler
    update_time!(get_base(cellbuffers), Δt)
    scatter!(cellbuffers)
    scatter!(assemblers)
    num_tasks = Threads.nthreads()
    for saved_set_chunks in buffer.saved_set_chunks
        set_chunks = ChunkIterator(saved_set_chunks)
        Base.Experimental.@sync begin 
            for taskid in 1:num_tasks
                cellbuffer = get_local(cellbuffers, taskid)
                assembler = get_local(assemblers, taskid)
                Threads.@spawn begin
                    for set_chunk in set_chunks
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

function threaded_assemble_domain4!(assemblers::TaskLocals, new_states, buffer::ThreadedDomainBuffer; a=nothing, aold=nothing, old_states=nothing, Δt=NaN)
    cellbuffers = buffer.cellbuffers #::TaskLocals
    sdh = buffer.sdh # SubDofHandler
    update_time!(get_base(cellbuffers), Δt)
    scatter!(cellbuffers)
    scatter!(assemblers)
    num_tasks = Threads.nthreads()
    for saved_set_chunks in buffer.saved_set_chunks
        set_chunks = ChunkIterator(saved_set_chunks)
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