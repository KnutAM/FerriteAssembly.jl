function work!(worker, buffer::Union{AbstractDomainBuffer, DomainBuffers}; a = nothing, aold = nothing)
    return work!(worker, Simulation(buffer, a, aold))
end

"""
    work!(worker, sim::Simulation, [coupled_simulations::CoupledSimulations])

Perform the work according to `worker` over the domain(s) in `sim`.

**Advance usage:** By passing the optional `coupled_simulations`, values from those simulations 
(e.g. state variables and local dof-values) become available on the local level via 
[`get_coupled_buffer`](@ref). This requires that the domainbuffer(s) in `sim` has been coupled 
using [`couple_buffers`](@ref).

    work!(worker, db::Union{AbstractDomainBuffer, Dict}; a = nothing, aold = nothing)

Simplified interface that doesn't support coupled simulations, directly forwarded to 
`work!(worker, Simulation(db, a, aold))`. The global degree of freedom vectors, `a` and `aold`,
make their corresponding local values available. If not passed, the local values are `NaN`s.
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
