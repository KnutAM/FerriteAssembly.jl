function work!(worker, buffer::Union{AbstractDomainBuffer, DomainBuffers}; a = nothing, aold = nothing)
    return work!(worker, Simulation(buffer, a, aold))
end

"""
    work!(worker, sim::Simulation, [coupled_simulations::CoupledSimulations])

Perform the work according to `worker` over the domain(s) in `sim`.

**Advance usage:** By passing the optional `coupled_simulations`, values from those simulations
(e.g. state variables and local dof-values) become available on the local level via
[`get_coupled_buffer`](@ref). The coupled buffer link is established fresh on each `work!` call
(no separate setup-time coupling step is needed), and always reflects whichever simulation is
currently supplied.

For threaded work, coupling reuses each partner's own per-task buffers directly (no partner
buffer content is copied) - so the partner's own task count must equal the primary domain's task
count (an `ArgumentError` is thrown otherwise); a sequential partner counts as having 1 task. A
small, fixed-size (task × partner count, not cell count) linking wrapper is still (re)constructed
each `work!` call.

!!! warning "Coupled simulations must not be mutated or shared concurrently"
    A coupled partner's degree-of-freedom vectors and state are read directly (not copied) while
    reiniting the per-cell/per-task views. This is race-free with respect to *this* `work!` call,
    but: (1) the coupled simulation(s) must not be assembled/updated by another, concurrently-
    running `work!` call while being read here - e.g. staggered solves must alternate (solve one
    domain, then the other), never run both domains' `work!` calls at the same time; and (2) for
    threaded work, the *same* coupled simulation must not be read by two different, concurrently-
    running `work!` calls either - since each reuses the partner's per-task buffers directly
    (primary task `i` writes into partner task `i`'s buffer), two concurrent primaries coupled to
    the same partner would race on those same buffers.

    work!(worker, db::Union{AbstractDomainBuffer, Dict}; a = nothing, aold = nothing)

Simplified interface that doesn't support coupled simulations, directly forwarded to
`work!(worker, Simulation(db, a, aold))`. The global degree of freedom vectors, `a` and `aold`,
make their corresponding local values available. If not passed, the local values are `NaN`s.
"""
function work!(worker, multisim::MultiDomainSim, coupled_simulations = CoupledSimulations())
    for (name, sim) in multisim
        skip_this_domain(worker, name) && continue
        coupled = get_domain_simulation(coupled_simulations, name)
        work_domain_sequential!(worker, sim, coupled)
    end
end
function work!(worker, sim::SingleDomainSim, coupled_simulations = CoupledSimulations())
    work_domain_sequential!(worker, sim, coupled_simulations)
end
function work!(worker, multisim::MultiDomainThreadedSim, coupled_simulations = CoupledSimulations())
    if can_thread(worker)
        workers = TaskLocals(worker, num_tasks = get_num_tasks(multisim))
        for (name, sim) in multisim
            skip_this_domain(worker, name) && continue
            coupled = get_domain_simulation(coupled_simulations, name)
            work_domain_threaded!(workers, sim, coupled)
        end
    else
        for (name, sim) in multisim
            skip_this_domain(worker, name) && continue
            coupled = get_domain_simulation(coupled_simulations, name)
            work_domain_sequential!(worker, sim, coupled)
        end
    end
end
function work!(worker, sim::SingleDomainThreadedSim, coupled_simulations = CoupledSimulations())
    if can_thread(worker)
        workers = TaskLocals(worker; num_tasks = get_num_tasks(sim))
        work_domain_threaded!(workers, sim, coupled_simulations)
    else
        work_domain_sequential!(worker, sim, coupled_simulations)
    end
end

function work_domain_sequential!(worker, sim::Simulation{<:AbstractDomainBuffer}, coupled)
    itembuffer = get_base(get_itembuffer(sim)) # get_base if threaded buffer
    coupled_itembuffers = map(get_base, (get_itembuffer(coupled))) # NamedTuple
    cb = couple_itembuffers(itembuffer, coupled_itembuffers)
    for itemnr in getset(sim)
        reinit_buffer!(cb, sim, coupled, itemnr)
        work_single!(worker, cb)
    end
end

function work_domain_threaded!(workers, sim::SingleDomainThreadedSim, coupled)
    itembuffers = get_itembuffer(sim) #::TaskLocals
    num_tasks = get_num_tasks(sim) # Default to Threads.nthreads()
    for coupled_sim in coupled.sims
        @assert get_num_tasks(coupled_sim) == num_tasks
    end
    coupled_itembuffers = get_itembuffer(coupled) #::NamedTuple{keys, <:TaskLocals}
    scatter!(itembuffers)
    scatter!(workers)
    scatter!.(values(coupled_itembuffers))
    for chunk_vector in get_chunks(sim)
        taskchunks = TaskChunks(chunk_vector)
        Base.Experimental.@sync begin
            for taskid in 1:num_tasks
                itembuffer = get_local(itembuffers, taskid)
                coupled_itembuffer = get_local(itembuffers, taskid)
                cib = couple_itembuffers(itembuffer, coupled_itembuffer)
                worker = get_local(workers, taskid)
                Threads.@spawn begin
                    while true
                        taskchunk = get_chunk(taskchunks) # Vector{Int}
                        isempty(taskchunk) && break
                        for itemnr in taskchunk
                            reinit_buffer!(cib, sim, coupled, itemnr) # also reinits the linked coupled buffers
                            work_single!(worker, cib)
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
