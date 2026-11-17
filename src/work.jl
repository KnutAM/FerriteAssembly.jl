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

    work!(worker, db::Union{AbstractDomainBuffer, FerriteAssembly.DomainBuffers}; a = nothing, aold = nothing)

Simplified interface, directly forwarded to `work!(worker, Simulation(db, a, aold))`.
The global degree of freedom vectors, `a` and `aold`, make their corresponding local values
available. If not passed, the local values are `NaN`s.

For a `FerriteAssembly.DomainBuffers`, all domains must be either all sequential
or all threaded; mixing the two in the same dictionary is not supported.

**Threaded error behavior:** if an element routine throws during threaded work, the other
spawned tasks are not cancelled, but all of them are finished (not left running in the
background) before the exception propagates out of `work!`.

**Algebraic variables:** A threaded domain with algebraic variables (see the `algebraic_variables`
keyword of [`DomainSpec`](@ref)) is only worked multithreaded if the `worker` supports this, see
[`FerriteAssembly.can_thread_shared_dofs`](@ref). Otherwise, that domain is worked sequentially.
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
            if can_thread_domain(worker, sim)
                work_domain_threaded!(workers, sim)
            else
                work_domain_sequential!(worker, sim)
            end
        end
    else
        for (name, sim) in multisim
            skip_this_domain(worker, name) && continue
            work_domain_sequential!(worker, sim)
        end
    end
end
function work!(worker, sim::AbstractSingleDomainThreadedSim)
    if can_thread(worker) && can_thread_domain(worker, sim)
        workers = TaskLocals(worker; num_tasks = get_num_tasks(sim))
        work_domain_threaded!(workers, sim)
    else
        work_domain_sequential!(worker, sim)
    end
end

# The dofs of algebraic variables are the same for all items in a domain. These are
# updated once, here, in the base itembuffer(s), before working the domain.
function refresh_algebraic_dofs!(sim::Simulation{<:AbstractDomainBuffer})
    update_algebraic_dofs!(get_base(get_itembuffer(sim)), get_dofhandler(sim))
end
function refresh_algebraic_dofs!(sim::CoupledSimulation{<:AbstractDomainBuffer})
    refresh_algebraic_dofs!(sim.sim)
    map(refresh_algebraic_dofs!, sim.partners)
    return nothing
end

has_algebraic_dofs(sim::AbstractSimulation{<:AbstractDomainBuffer}) = has_algebraic_dofs(get_base(get_itembuffer(sim)))

# Can the `worker`, which fulfills `can_thread(worker)`, work the domain in `sim` multithreaded?
function can_thread_domain(worker, sim::AbstractSimulation{<:AbstractDomainBuffer})
    has_algebraic_dofs(sim) || return true
    can_thread_shared_dofs(worker) && return true
    @warn(string("The domain has algebraic variables, whose degrees of freedom are shared between all items. ",
        "The worker, ", nameof(typeof(worker)), ", does not support multithreaded work for such domains ",
        "(see `FerriteAssembly.can_thread_shared_dofs`), and the domain is worked sequentially. ",
        "For assembly, use `KeReAssembler(K, r; atomic = true)` to enable multithreading."),
        maxlog = 1, _id = Symbol(:shared_dofs_sequential_, nameof(typeof(worker))))
    return false
end

function work_domain_sequential!(worker, sim::AbstractSimulation{<:AbstractDomainBuffer})
    refresh_algebraic_dofs!(sim)
    itembuffer = get_base(get_itembuffer(sim)) # get_base if threaded buffer
    for itemnr in get_set(sim)
        reinit_buffer!(itembuffer, sim, itemnr)
        work_single!(worker, itembuffer)
    end
end

function work_domain_threaded!(workers, sim::AbstractSingleDomainThreadedSim)
    refresh_algebraic_dofs!(sim)
    scatter!(sim) # Includes scatter of the `itembuffers`
    itembuffers = get_itembuffer(sim) #::TaskLocals
    scatter!(workers)
    num_tasks = get_num_tasks(sim) # Default to Threads.nthreads()
    for chunk_vector in get_chunks(sim)
        taskchunks = TaskChunks(chunk_vector)
        spawn_tasks_safe_wait(num_tasks) do taskid, abort
            itembuffer = get_local(itembuffers, taskid)
            worker = get_local(workers, taskid)
            while !abort[]
                taskchunk = get_chunk(taskchunks) # Union{Vector{Int}, Nothing}
                taskchunk === nothing && break
                for itemnr in taskchunk
                    abort[] && break
                    reinit_buffer!(itembuffer, sim, itemnr)
                    work_single!(worker, itembuffer)
                end # itemnr
            end #chunk
        end
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
    can_thread_shared_dofs(worker)::Bool

Does the `worker` support multithreaded work on a domain with degrees of freedom that are shared
between all items in the domain? This is the case for domains with algebraic variables
(see the `algebraic_variables` keyword of [`DomainSpec`](@ref)). Coloring does not prevent race
conditions for such degrees of freedom. Defaults to `false`. Only used if `can_thread(worker)`,
and if `false`, domains with algebraic variables are worked sequentially.

A worker that only accumulates into task local values (or values that are unique for each item)
can return `true`. A worker that writes to global values indexed by the degrees of freedom
(e.g. assemblers) can only return `true` if these writes are safe to do concurrently,
e.g. by using atomic operations.
"""
can_thread_shared_dofs(::Any) = false

"""
    skip_this_domain(worker, name::String)

Should the domain with key `name` be skipped during work? Defaults to `false`.
Can be used to e.g. only loop over parts of a domain.
"""
skip_this_domain(::Any, ::String) = false # opt-in to skip domains (used for integration)
