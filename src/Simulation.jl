# Simulation 
struct Simulation{DB <: Union{DomainBuffers, AbstractDomainBuffer}, TA, TAO}
    db::DB
    a::TA
    aold::TAO
end

const SingleDomainSim = Simulation{<:DomainBuffer}
const MultiDomainSim = Simulation{<:Dict{String, <:DomainBuffer}}
const SingleDomainThreadedSim = Simulation{<:ThreadedDomainBuffer}
const MultiDomainThreadedSim = Simulation{<:Dict{String, <:ThreadedDomainBuffer}}

# Iterator interface
@inline function _iterate(sim::Simulation{<:Dict}, iter)
    iter === nothing && return nothing
    ((name, db), state) = iter              
    return ((name, Simulation(db, sim.a, sim.aold)), state)
end
Base.iterate(sim::Simulation{<:Dict}) = _iterate(sim, iterate(sim.db))
Base.iterate(sim::Simulation{<:Dict}, iter) = _iterate(sim, iterate(sim.db, iter))

get_itembuffer(sim::Simulation{<:AbstractDomainBuffer}) = get_itembuffer(sim.db)
get_num_tasks(sim::Simulation{<:AbstractDomainBuffer}) = get_num_tasks(sim.db)
get_chunks(sim::Simulation{<:AbstractDomainBuffer}) = get_chunks(sim.db)
getset(sim::Simulation{<:AbstractDomainBuffer}) = getset(sim.db)

get_dofhandler(sim::Simulation) = get_dofhandler(sim.db)
get_state(sim::Simulation, args::Vararg{Any,N}) where N = get_state(sim.db, args...)
get_old_state(sim::Simulation, args::Vararg{Any,N}) where N = get_old_state(sim.db, args...)

get_domain(sim::Simulation{<:Dict}, name::String) = Simulation(sim.db[name], sim.a, sim.aold)

struct CoupledSimulations{NT <: NamedTuple{<:Any, <:NTuple{<:Any, Simulation}}}
    sims::NT
end
CoupledSimulations(; kwargs...) = CoupledSimulations(NamedTuple{keys(kwargs)}(values(kwargs)))

function get_domain(cs::CoupledSimulations, name::String)
    return CoupledSimulations(map(s -> get_domain(s, name), cs.sims))
end
