function work!(worker, buffer; a = nothing, aold = nothing)
    return work!(worker, Simulation(buffer, a, aold))
end

"""
    work!(worker, simulation, [coupled_simulations])

Perform the work according to `worker` over the domain(s) specified by 
`buffer`. Current, `a`, and old, `aold`, global degree of freedom vectors 
are passed to get these values passed into the innermost user-defined functions.
If not passed (or as `nothing`), `NaN` values are passed into the innermost functions. 
"""
function work!(worker, multisim::MultiDomainSim, coupled_simulations = CoupledSimulations())
    for (name, sim) in multisim
        skip_this_domain(worker, name) && continue
        coupled = get_domain(coupled_simulations, name)
        work_domain_sequential!(worker, sim, coupled)
    end
end
function work!(worker, sim::SingleDomainSim, coupled_simulations = CoupledSimulations())
    work_domain_sequential!(worker, sim, coupled_simulations)
end
function work!(worker, multisim::MultiDomainThreadedSim, coupled_simulations = CoupledSimulations())
    if can_thread(worker)
        workers = TaskLocals(worker)
        for (name, sim) in multisim
            skip_this_domain(worker, name) && continue
            coupled = get_domain(coupled_simulations, name)
            work_domain_threaded!(workers, sim, coupled)
        end
    else
        for (name, sim) in multisim
            skip_this_domain(worker, name) && continue
            coupled = get_domain(coupled_simulations, name)
            work_domain_sequential!(worker, sim, coupled)
        end
    end
end
function work!(worker, sim::SingleDomainThreadedSim, coupled_simulations = CoupledSimulations())
    if can_thread(worker)
        workers = TaskLocals(worker)
        work_domain_threaded!(workers, sim, coupled_simulations)
    else
        work_domain_sequential!(worker, sim, coupled_simulations)
    end
end

function work_domain_sequential!(worker, sim::Simulation{<:AbstractDomainBuffer}, coupled)
    itembuffer = get_base(get_itembuffer(sim)) # get_base if threaded buffer
    for itemnr in getset(sim)
        reinit_buffer!(itembuffer, sim, coupled, itemnr)
        work_single!(worker, itembuffer)
    end
end

function work_domain_threaded!(workers, sim::SingleDomainThreadedSim, coupled)
    itembuffers = get_itembuffer(sim) #::TaskLocals
    scatter!(itembuffers)
    scatter!(workers)
    num_tasks = get_num_tasks(sim) # Default to Threads.nthreads()
    for chunk_vector in get_chunks(sim)
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
                            reinit_buffer!(itembuffer, sim, coupled, itemnr)
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
