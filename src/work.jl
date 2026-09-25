function work!(worker, buffer::Union{AbstractDomainBuffer, DomainBuffers}; a = nothing, aold = nothing)
    return work!(worker, Simulation(buffer, a, aold))
end

"""
    work!(worker, sim::Simulation)

Perform the work according to `worker` over the domain(s) in `sim`.

**Coupled simulations:** To make values from other simulations (e.g. state variables and
local dof-values) available on the local level via [`get_coupled_buffer`](@ref), build a
[`CoupledSimulations`](@ref) group and call `work!(worker, group.member_name)` instead;
the member handle already carries its resolved coupling.

    work!(worker, db::Union{AbstractDomainBuffer, Dict}; a = nothing, aold = nothing)

Simplified interface, directly forwarded to `work!(worker, Simulation(db, a, aold))`.
The global degree of freedom vectors, `a` and `aold`, make their corresponding local values
available. If not passed, the local values are `NaN`s.
"""
function work!(worker, multisim::AbstractMultiDomainSim)
    for (name, sim) in multisim
        skip_this_domain(worker, name) && continue
        work_domain_sequential!(worker, sim)
    end
end
function work!(worker, sim::AbstractSingleDomainSim)
    work_domain_sequential!(worker, sim)
end
function work!(worker, multisim::AbstractMultiDomainThreadedSim)
    isempty(get_domainbuffer(multisim)) && return nothing
    if can_thread(worker)
        workers = TaskLocals(worker, num_tasks = get_num_tasks(multisim))
        for (name, sim) in multisim
            skip_this_domain(worker, name) && continue
            work_domain_threaded!(workers, sim)
        end
    else
        for (name, sim) in multisim
            skip_this_domain(worker, name) && continue
            work_domain_sequential!(worker, sim)
        end
    end
end
function work!(worker, sim::AbstractSingleDomainThreadedSim)
    if can_thread(worker)
        workers = TaskLocals(worker; num_tasks = get_num_tasks(sim))
        work_domain_threaded!(workers, sim)
    else
        work_domain_sequential!(worker, sim)
    end
end

function work_domain_sequential!(worker, sim::AbstractSimulation{<:AbstractDomainBuffer})
    itembuffer = get_base(get_itembuffer(sim)) # get_base if threaded buffer
    for itemnr in getset(sim)
        reinit_buffer!(itembuffer, sim, itemnr)
        work_single!(worker, itembuffer)
    end
end

function work_domain_threaded!(workers, sim::AbstractSingleDomainThreadedSim)
    scatter!(sim) # Includes scatter of the `itembuffers`
    itembuffers = get_itembuffer(sim) #::TaskLocals
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
                        taskchunk = get_chunk(taskchunks) # Union{Vector{Int}, Nothing}
                        taskchunk === nothing && break
                        for itemnr in taskchunk
                            reinit_buffer!(itembuffer, sim, itemnr)
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
