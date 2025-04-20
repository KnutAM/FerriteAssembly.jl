"""
    work!(worker, buffer; a=nothing, aold=nothing)

Perform the work according to `worker` over the domain(s) specified by 
`buffer`. Current, `a`, and old, `aold`, global degree of freedom vectors 
are passed to get these values passed into the innermost user-defined functions.
If not passed (or as `nothing`), `NaN` values are passed into the innermost functions. 
"""
function work!(worker, buffers::Dict{String,<:DomainBuffer}; kwargs...)
    for (name, buffer) in buffers
        skip_this_domain(worker, name) && continue
        work_domain_sequential!(worker, buffer; kwargs...)
    end
end
function work!(worker, buffer::DomainBuffer; kwargs...)
    work_domain_sequential!(worker, buffer; kwargs...)
end
function work!(worker, buffers::Dict{String,<:ThreadedDomainBuffer}; kwargs...)
    if can_thread(worker)
        workers = TaskLocals(worker)
        for (name, buffer) in buffers
            skip_this_domain(worker, name) && continue
            work_domain_threaded!(workers, buffer; kwargs...)
        end
    else
        for (name, buffer) in buffers
            skip_this_domain(worker, name) && continue
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
    itembuffer = get_base(get_itembuffer(domainbuffer)) # get_base if threaded buffer
    for itemnr in getset(domainbuffer)
        reinit_buffer!(itembuffer, domainbuffer, itemnr; kwargs...)
        work_single!(worker, itembuffer)
    end
end

function work_domain_threaded!(workers, domainbuffer; kwargs...)
    itembuffers = get_itembuffer(domainbuffer) #::TaskLocals
    scatter!(itembuffers)
    scatter!(workers)
    num_tasks = get_num_tasks(domainbuffer) # Default to Threads.nthreads()
    for chunk_vector in get_chunks(domainbuffer)
        taskchunks = TaskChunks(chunk_vector)
        Base.Experimental.@sync begin 
            for taskid in 1:num_tasks
                itembuffer = get_local(itembuffers, taskid)
                worker = get_local(workers, taskid)
                Threads.@spawn begin
                    while true
                        taskchunk = get_chunk(taskchunks) # Vector{Int}
                        isempty(taskchunk) && break
                        for itemnr in taskchunk
                            reinit_buffer!(itembuffer, domainbuffer, itemnr; kwargs...)
                            work_single!(worker, itembuffer)
                        end # itemnr
                    end #chunk
                end #spawn
            end #taskid
        end #sync
    end #chunk_vectors
    gather!(itembuffers)
    gather!(workers)
end

# Worker interface
"""
    can_thread(worker)::Bool

Does the worker support multithreaded work? Defaults to `false`. 
If this returns `true`, the worker must support the `TaskLocals` interface. 
"""
can_thread(::Any) = false

"""
    skip_this_domain(worker, name::String)

Should the domain with key `name` be skipped during work? Defaults to `false`.
Can be used to e.g. only loop over parts of a domain.
"""
skip_this_domain(::Any, ::String) = false # opt-in to skip domains (used for integration)
