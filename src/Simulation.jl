abstract type AbstractSimulation{DB} end

const AbstractSingleDomainSim = AbstractSimulation{<:DomainBuffer}
const AbstractMultiDomainSim = AbstractSimulation{<:Dict{String, <:DomainBuffer}}
const AbstractSingleDomainThreadedSim = AbstractSimulation{<:ThreadedDomainBuffer}
const AbstractMultiDomainThreadedSim = AbstractSimulation{<:Dict{String, <:ThreadedDomainBuffer}}

# Must be defined
"""
    get_domainbuffer(sim::AbstractSimulation)

Accessor for the automatic forwarding for domainbuffer methods to work
"""
function get_domainbuffer end

# Forwarding for public API
get_material(sim::AbstractSimulation, args::Vararg{Any, N}) where N = get_material(get_domainbuffer(sim), args...)
get_dofhandler(sim::AbstractSimulation) = get_dofhandler(get_domainbuffer(sim))
get_grid(sim::AbstractSimulation) = get_grid(get_domainbuffer(sim))
get_state(sim::AbstractSimulation, args::Vararg{Any, N}) where N = get_state(get_domainbuffer(sim), args...)
get_old_state(sim::AbstractSimulation, args::Vararg{Any, N}) where N = get_old_state(get_domainbuffer(sim), args...)
getset(sim::AbstractSimulation, args::Vararg{Any, N}) where N = getset(get_domainbuffer(sim), args...)
update_states!(sim::AbstractSimulation; kwargs...) = update_states!(get_domainbuffer(sim); kwargs...)
set_time_increment!(sim::AbstractSimulation, Δt) = set_time_increment!(get_domainbuffer(sim), Δt)
revert_states!(sim::AbstractSimulation) = revert_states!(get_domainbuffer(sim))

# Forwarding for internal API
get_num_tasks(sim::AbstractSimulation) = get_num_tasks(get_domainbuffer(sim))
get_chunks(sim::AbstractSimulation{<:AbstractDomainBuffer}) = get_chunks(get_domainbuffer(sim))
get_itembuffer(sim::AbstractSimulation, args::Vararg{Any, N}) where {N} = get_itembuffer(get_domainbuffer(sim), args...)


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
        } <: AbstractSimulation{DB}
    db::DB
    a::TA
    aold::TAO
end
Simulation(db::Union{DomainBuffers, AbstractDomainBuffer}, a = nothing, aold = nothing) = Simulation(db, a, aold)

get_domainbuffer(sim::Simulation) = sim.db

## Iterator interface
@inline function _iterate(sim::Simulation{<:DomainBuffers}, iter)
    iter === nothing && return nothing
    ((name, db), state) = iter
    return ((name, Simulation(db, sim.a, sim.aold)), state)
end
Base.iterate(sim::Simulation{<:DomainBuffers}) = _iterate(sim, iterate(sim.db))
Base.iterate(sim::Simulation{<:DomainBuffers}, iter) = _iterate(sim, iterate(sim.db, iter))

"""
    CoupledSimulation(sim, partners)

A handle to one primary member of a [`CoupledSimulations`](@ref) group (e.g. `group.a`).
`sim` is a [`Simulation`](@ref) whose domain buffer(s) have been rebuilt with coupled
itembuffers. `partners` holds the resolved partner `Simulation`s this member reads from: a
`NamedTuple{name}` of partner `Simulation`s for a single-domain member, or a
`Dict{String,<:NamedTuple}` (one `NamedTuple` of partner `Simulation`s per domain name) for a
multi-domain member. This is the single, canonical copy of that information — passed into
[`reinit_buffer!`](@ref) at call time rather than duplicated into every task-local buffer.

Forwards the ordinary [`Simulation`](@ref) accessor API (`.a`, `.aold`, `.db`,
`get_dofhandler`, `get_state`, `set_time_increment!`, `update_states!`, etc.).
"""
struct CoupledSimulation{DB, S <: Simulation{DB}, P} <: AbstractSimulation{DB}
    sim::S
    partners::P
end

get_domainbuffer(sim::CoupledSimulation) = get_domainbuffer(sim.sim)

replace_material(::CoupledSimulation, args...; kwargs...) = throw(ArgumentError(
    "replace_material on a CoupledSimulations member is not supported; use " *
    "replace_material(group, member_name, f) to rebuild the whole group instead."))

# Per-domain iteration for a multi-domain member, mirroring `Simulation{<:DomainBuffers}`'s
# own iteration but pairing each per-domain `Simulation` with its own slice of `partners`.
function Base.iterate(csim::CoupledSimulation{<:DomainBuffers})
    it = iterate(csim.sim)
    it === nothing && return nothing
    ((name, dsim), st) = it
    return ((name, CoupledSimulation(dsim, csim.partners[name])), st)
end
function Base.iterate(csim::CoupledSimulation{<:DomainBuffers}, st)
    it = iterate(csim.sim, st)
    it === nothing && return nothing
    ((name, dsim), st2) = it
    return ((name, CoupledSimulation(dsim, csim.partners[name])), st2)
end

_scatter_partner_container!(c::TaskLocals) = scatter!(c)
_scatter_partner_container!(::Any) = nothing

_flatten_partner_sims(partners::NamedTuple) = values(partners)
_flatten_partner_sims(partners_by_domain::Dict) = (psim for nt in values(partners_by_domain) for psim in values(nt))

# Scatter every reachable partner's task-local buffers from its base once, before any
# per-cell work, so a threaded reader always observes the partner's *current* state (e.g. its
# time increment) even if the partner itself has not been `work!`ed since it last changed.
function _scatter_all_partners!(csim::CoupledSimulation)
    for psim in _flatten_partner_sims(csim.partners)
        _scatter_partner_container!(get_itembuffer(psim.db))
    end
    return nothing
end

# Hook hit once at the start of every top-level `work!` call (see `work.jl`); a no-op for a
# plain `Simulation`, overridden here so a `CoupledSimulation`'s partners are scattered before
# any per-cell work, without `work.jl` needing to know coupling exists.
_prepare_work!(csim::CoupledSimulation) = _scatter_all_partners!(csim)
