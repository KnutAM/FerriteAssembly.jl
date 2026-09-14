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

