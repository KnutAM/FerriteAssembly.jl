"""
    Simulation(db, a = nothing, aold = nothing)

A `Simulation` is a collection of the simulation domain(s) `db`, and the 
global degree of freedom vectors, `a` and `aold`. 

**Note:** 
If `a` or `aold` are not provided, the local vectors will have `NaN` values.
"""
struct Simulation{
        DB  <: Union{DomainBuffers, AbstractDomainBuffer}, 
        TA  <: Union{Nothing, AbstractVector}, 
        TAO <: Union{Nothing, AbstractVector}
        }
    db::DB
    a::TA
    aold::TAO
end
Simulation(db::Union{DomainBuffers, AbstractDomainBuffer}, a = nothing, aold = nothing) = Simulation(db, a, aold)

const SingleDomainSim = Simulation{<:DomainBuffer}
const MultiDomainSim = Simulation{<:Dict{String, <:DomainBuffer}}
const SingleDomainThreadedSim = Simulation{<:ThreadedDomainBuffer}
const MultiDomainThreadedSim = Simulation{<:Dict{String, <:ThreadedDomainBuffer}}

# Forwarding for public API
get_material(sim::Simulation, args::Vararg{Any, N}) where N = get_material(sim.db, args...)
get_dofhandler(sim::Simulation) = get_dofhandler(sim.db)
get_grid(sim::Simulation) = get_grid(sim.db)
get_state(sim::Simulation, args::Vararg{Any, N}) where N = get_state(sim.db, args...)
get_old_state(sim::Simulation, args::Vararg{Any, N}) where N = get_old_state(sim.db, args...)
getset(sim::Simulation, args::Vararg{Any, N}) where N = getset(sim.db, args...)
update_states!(sim::Simulation; kwargs...) = update_states!(sim.db; kwargs...)
set_time_increment!(sim::Simulation, Δt) = set_time_increment!(sim.db, Δt)
revert_states!(sim::Simulation) = revert_states!(sim.db)

# Forwarding for internal API
get_num_tasks(sim::Simulation) = get_num_tasks(sim.db)
get_chunks(sim::Simulation{<:AbstractDomainBuffer}) = get_chunks(sim.db)
get_itembuffer(sim::Simulation, args::Vararg{Any, N}) where {N} = get_itembuffer(sim.db, args...)

# Internal API
get_domain_simulation(sim::Simulation{<:DomainBuffers}, name::String) = Simulation(sim.db[name], sim.a, sim.aold)
## Iterator interface
@inline function _iterate(sim::Simulation{<:DomainBuffers}, iter)
    iter === nothing && return nothing
    ((name, db), state) = iter              
    return ((name, Simulation(db, sim.a, sim.aold)), state)
end
Base.iterate(sim::Simulation{<:DomainBuffers}) = _iterate(sim, iterate(sim.db))
Base.iterate(sim::Simulation{<:DomainBuffers}, iter) = _iterate(sim, iterate(sim.db, iter))

"""
    CoupledSimulations(; key1 = sim1::Simulation, key2 = sim2::Simulation, ...)

Setup the collection of coupled simulations to allow values (such as state variables and
local dof-values) from these simulations to be available when `work!`ing another simulation.
The coupled itembuffer on the local level is accessed with [`get_coupled_buffer`](@ref).
"""
struct CoupledSimulations{NT <: NamedTuple{<:Any, <:NTuple{<:Any, Simulation}}}
    sims::NT
end
CoupledSimulations(; kwargs...) = CoupledSimulations(NamedTuple{keys(kwargs)}(values(kwargs)))

"""
    get_itembuffer(coupled::CoupledSimulations, num_tasks::Int)

Return a `NamedTuple` (same keys as `coupled.sims`) of each coupled partner's own itembuffer
(its `TaskLocals`, if threaded, or its single buffer, if sequential) - reused directly, so no
partner buffer content is copied. Each partner's own task count must equal `num_tasks`, the
primary (coupling) domain's own task count: reusing a partner's per-task buffers is only
race-free if primary task `i` always maps to partner task `i`, one-to-one, for every task - which
also means the same partner must not be coupled to by two different, concurrently-running `work!`
calls (see the warning on [`work!`](@ref)).
"""
function get_itembuffer(coupled::CoupledSimulations, num_tasks::Int)
    ks = keys(coupled.sims)
    return NamedTuple{ks}(map(ks) do k
        sim = coupled.sims[k]
        partner_tasks = get_num_tasks(sim)
        partner_tasks == num_tasks || throw(ArgumentError(
            "Coupled simulation `:$k` has $partner_tasks task(s), but the primary domain has " *
            "$num_tasks. These must match for threaded coupled work: set matching `num_tasks` " *
            "when setting up both domains (a sequential domain counts as 1 task)."
        ))
        get_itembuffer(sim)
    end)
end

function get_domain_simulation(cs::CoupledSimulations, name::String)
    # Need to return a named tuple with only the simulations that have a domain called `name`
    sims = Pair{Symbol, Simulation}[]
    for (key, sim) in zip(keys(cs.sims), values(cs.sims))
        if haskey(sim.db, name)
            push!(sims, key => get_domain_simulation(sim, name))
        end
    end
    return CoupledSimulations(NamedTuple(sims))
#    return CoupledSimulations(map(s -> get_domain_simulation(s, name), cs.sims))
end