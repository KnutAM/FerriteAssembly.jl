# CoupledSimulations group construction and the cell-buffer-dependent parts of coupling.
#
# CoupledSimulation itself (struct, accessors, iteration, partner-scatter hook) lives in
# Simulation.jl since it needs no cell-buffer types. Everything here references CellBuffer,
# CoupledCellBuffer, or AutoDiffCellBuffer, so this file must be included after those exist
# (see FerriteAssembly.jl's include order).

"""
    reinit_buffer!(cb::CoupledCellBuffer, sim::CoupledSimulation, cellnum::Int)

Reinitialize the reader's own `cb.primary` against `sim.sim`, then reinitialize each partner
buffer in `cb.partner_buffers` against its own partner `Simulation` stored in `sim.partners`.
Partner reinitialization does not recurse: partner buffers are plain `CellBuffer`s, so no
further coupling initialization happens.
"""
function reinit_buffer!(cb::CoupledCellBuffer, sim::CoupledSimulation, cellnum::Int)
    reinit_buffer!(cb.primary, getfield(sim, :sim), cellnum)
    map((b, s) -> (reinit_buffer!(b, s, cellnum); nothing), cb.partner_buffers, sim.partners)
    return nothing
end

struct CoupledSimulations{P<:NamedTuple, R<:NamedTuple, M<:NamedTuple}
    primaries::P
    refs::R
    members::M
end

"""
    CoupledSimulations(primaries::NamedTuple; refs::NamedTuple = NamedTuple())

Build a group of mutually-wired simulations from `primaries` (members that read partners
and are worked via the group) and, optionally, `refs` (members with no outgoing
dependencies, still accessible/workable through the group but never rewired themselves).

Each primary reads every other primary and every ref (excluded: itself). Names must be
unique across `primaries` and `refs`. Member access is direct/nonrecursive: `g.a`'s view of
`g.b` exposes `b`'s own local values, not `b`'s further coupling.

```julia
g = CoupledSimulations((a = sima, b = simb, c = simc))              # mutual
g = CoupledSimulations((a = sima,); refs = (b = simb,))             # one-way: a reads b
g = CoupledSimulations((a = sima, b = simb); refs = (c = simc,))    # mixed

work!(worker_a, g.a)
```

See the package documentation for the full setup-validation and replacement contract.
"""
function CoupledSimulations(primaries::NamedTuple; refs::NamedTuple = NamedTuple())
    isempty(primaries) && throw(ArgumentError("`primaries` must be a nonempty named tuple of `Simulation`s"))
    all(v -> v isa Simulation, primaries) || throw(ArgumentError("`primaries` values must be `Simulation`s"))
    all(v -> v isa Simulation, refs) || throw(ArgumentError("`refs` values must be `Simulation`s"))
    overlap = intersect(keys(primaries), keys(refs))
    isempty(overlap) || throw(ArgumentError("primary and ref names must be unique, overlap: $overlap"))
    reserved = intersect(union(keys(primaries), keys(refs)), (:primaries, :refs, :members))
    isempty(reserved) || throw(ArgumentError(
        "member name(s) $reserved are reserved and would shadow `CoupledSimulations` internals"))

    all_members = merge(primaries, refs)
    validate_storage_identity(all_members)
    validate_task_counts_positive(all_members)
    members = NamedTuple{keys(primaries)}(
        Tuple(build_coupled_simulation(name, sim, all_members) for (name, sim) in pairs(primaries))
    )
    return CoupledSimulations(primaries, refs, members)
end

function Base.getproperty(cs::CoupledSimulations, name::Symbol)
    name in (:primaries, :refs, :members) && return getfield(cs, name)
    members = getfield(cs, :members)
    haskey(members, name) && return members[name]
    refs = getfield(cs, :refs)
    haskey(refs, name) && return refs[name]
    throw(ArgumentError("CoupledSimulations has no member named `$name`"))
end

Base.propertynames(cs::CoupledSimulations) = (:primaries, :refs, :members, keys(getfield(cs, :primaries))..., keys(getfield(cs, :refs))...)

unwrap_cb(cb::CellBuffer) = cb
unwrap_cb(ad::AutoDiffCellBuffer) = ad.cb

select_partner(p::TaskLocals, i::Int) = get_local(p, i)
select_partner(p, ::Int) = p

_is_autodiff(ib::AutoDiffCellBuffer) = true
_is_autodiff(ib) = ib isa TaskLocals && get_base(ib) isa AutoDiffCellBuffer

function build_coupled_itembuffer(reader_ibuf, partner_containers::NamedTuple)
    autodiff = _is_autodiff(reader_ibuf)
    wrap(primary_cb, partners_nt) = autodiff ?
        AutoDiffCellBuffer(CoupledCellBuffer(primary_cb, partners_nt)) :
        CoupledCellBuffer(primary_cb, partners_nt)
    if reader_ibuf isa TaskLocals
        n = length(get_locals(reader_ibuf))
        base = wrap(unwrap_cb(get_base(reader_ibuf)), map(unwrap_cb ∘ get_base, partner_containers))
        locals = [wrap(unwrap_cb(get_local(reader_ibuf, i)),
                        map(c -> unwrap_cb(select_partner(c, i)), partner_containers)) for i in 1:n]
        return TaskLocals(base, locals)
    else
        return wrap(unwrap_cb(reader_ibuf), map(unwrap_cb ∘ get_base, partner_containers))
    end
end

function validate_domain_pair(reader_db::AbstractDomainBuffer, partner_db::AbstractDomainBuffer, partner_name::Symbol)
    get_grid(reader_db) === get_grid(partner_db) || throw(ArgumentError(
        "coupling partner `$partner_name` uses a different grid than the reader"))
    for (role, db) in ((:reader, reader_db), (Symbol(partner_name), partner_db))
        ib = get_base(get_itembuffer(db))
        (ib isa CellBuffer || ib isa AutoDiffCellBuffer) || throw(ArgumentError(
            "coupling only supports `CellBuffer`/autodiff cell buffers, got $(typeof(ib)) for `$role`"))
    end
    issubset(getset(reader_db), getset(partner_db)) || throw(ArgumentError(
        "coupling partner `$partner_name` does not cover all cells read by the reader"))
    reader_threaded = reader_db isa ThreadedDomainBuffer
    if reader_threaded
        reader_tasks = get_num_tasks(reader_db)
        reader_tasks > 0 || throw(ArgumentError("task count must be positive"))
        partner_tasks = partner_db isa ThreadedDomainBuffer ? get_num_tasks(partner_db) : 1
        reader_tasks == partner_tasks || throw(ArgumentError(
            "threaded reader with $reader_tasks tasks requires coupling partner `$partner_name` to provide " *
            "$reader_tasks task-local buffers (a sequential partner counts as 1 slot); got $partner_tasks"))
    end
    return nothing
end

# Returns (new_db, partner_sims::NamedTuple): the rebuilt domain buffer with coupled
# itembuffer(s), and the resolved per-domain partner `Simulation`s (for the caller to store
# on the owning `CoupledSimulation`, not on the itembuffer itself).
function build_coupled_domain(reader_db::AbstractDomainBuffer, partners::NamedTuple)
    for (pname, p) in pairs(partners)
        validate_domain_pair(reader_db, p.db, pname)
    end
    reader_ibuf = get_itembuffer(reader_db)
    partner_containers = map(p -> get_itembuffer(p.db), partners)
    coupled_ibuf = build_coupled_itembuffer(reader_ibuf, partner_containers)
    new_db = setproperties(reader_db; itembuffer = coupled_ibuf)
    partner_sims = map(p -> p.sim, partners)
    return new_db, partner_sims
end

# Resolve, for a single reader domain (named `dname` when the reader is a `Dict`, or
# `nothing` for a single-domain reader), the single-domain `Simulation` of a partner
# (sharing the partner's own `a`/`aold`). Errors if a partner does not provide a required
# domain, or if reader/partner shapes are mixed.
function partner_domain_sim(dname::Union{Nothing,String}, partner_name::Symbol, partner_sim::Simulation)
    pdb = partner_sim.db
    if dname === nothing
        pdb isa DomainBuffers && throw(ArgumentError(
            "mixed single-domain/dictionary coupling is not supported (reader is single-domain, " *
            "partner `$partner_name` is a domain dictionary)"))
        return partner_sim
    else
        pdb isa DomainBuffers || throw(ArgumentError(
            "mixed single-domain/dictionary coupling is not supported (reader is a domain dictionary, " *
            "partner `$partner_name` is single-domain)"))
        haskey(pdb, dname) || throw(ArgumentError(
            "coupling partner `$partner_name` does not supply required domain \"$dname\""))
        return Simulation(pdb[dname], partner_sim.a, partner_sim.aold)
    end
end

function build_coupled_simulation(name::Symbol, reader_sim::Simulation, all_members::NamedTuple)
    partner_names = Tuple(k for k in keys(all_members) if k != name)
    partner_sims = NamedTuple{partner_names}(Tuple(all_members[k] for k in partner_names))
    reader_db = reader_sim.db
    if reader_db isa DomainBuffers
        if isempty(reader_db)
            new_db = reader_db # nothing to couple; preserves the original (correctly-typed) empty Dict
            partners_by_domain = Dict{String, NamedTuple}()
        else
            built = Any[]
            partners_by_domain = Dict{String, Any}()
            for (dname, rdb) in reader_db
                partners = NamedTuple{partner_names}(Tuple(
                    let psim_dom = partner_domain_sim(dname, pname, psim)
                        (sim = psim_dom, db = psim_dom.db)
                    end for (pname, psim) in pairs(partner_sims)
                ))
                ndb, dpartner_sims = build_coupled_domain(rdb, partners)
                push!(built, dname => ndb)
                partners_by_domain[dname] = dpartner_sims
            end
            new_db = Dict(built...) # infers the narrowest common concrete value type, matching MultiDomain(Threaded)Sim dispatch
        end
        new_sim = Simulation(new_db, reader_sim.a, reader_sim.aold)
        return CoupledSimulation(new_sim, partners_by_domain)
    else
        partners = NamedTuple{partner_names}(Tuple(
            let psim_dom = partner_domain_sim(nothing, pname, psim)
                (sim = psim_dom, db = psim_dom.db)
            end for (pname, psim) in pairs(partner_sims)
        ))
        new_db, dpartner_sims = build_coupled_domain(reader_db, partners)
        new_sim = Simulation(new_db, reader_sim.a, reader_sim.aold)
        return CoupledSimulation(new_sim, dpartner_sims)
    end
end

_scratch_identity(cb::CellBuffer) = cb.ae # survives replace_material's setproperties (fields copied by reference)

function _domain_entries(name::Symbol, sim::Simulation)
    db = sim.db
    db isa DomainBuffers && return [(name, dname, _scratch_identity(unwrap_cb(get_base(get_itembuffer(d))))) for (dname, d) in db]
    return [(name, "", _scratch_identity(unwrap_cb(get_base(get_itembuffer(db)))))]
end

function validate_storage_identity(all_members::NamedTuple)
    entries = reduce(vcat, (_domain_entries(name, sim) for (name, sim) in pairs(all_members)))
    for i in eachindex(entries), j in (i+1):length(entries)
        if entries[i][3] === entries[j][3]
            throw(ArgumentError(
                "members `$(entries[i][1])` (domain \"$(entries[i][2])\") and `$(entries[j][1])` " *
                "(domain \"$(entries[j][2])\") alias the same underlying item-buffer storage; " *
                "each member must own distinct mutable scratch"))
        end
    end
    return nothing
end

# Every member's own domain(s) must have a positive task count, independent of whether that
# member is ever paired as a "reader" against a partner (e.g. a ref, or a sole primary with
# no partners, would otherwise never be checked).
function _validate_own_task_count(name::Symbol, db::AbstractDomainBuffer)
    db isa ThreadedDomainBuffer || return nothing
    get_num_tasks(db) > 0 || throw(ArgumentError("member `$name` has a nonpositive task count"))
    return nothing
end

function validate_task_counts_positive(all_members::NamedTuple)
    for (name, sim) in pairs(all_members)
        db = sim.db
        if db isa DomainBuffers
            for (_, d) in db
                _validate_own_task_count(name, d)
            end
        else
            _validate_own_task_count(name, db)
        end
    end
    return nothing
end

"""
    replace_material(g::CoupledSimulations, member::Symbol, f; domain = nothing)

Return a new `CoupledSimulations` group in which `member`'s material has been replaced by
`f` (applied as `f(old_material)`), either for the whole member (`domain = nothing`) or only
for the named domain of a multi-domain member. Rebuilds the whole group (rerunning
constructor validation and autodiff configuration construction); other members are reused by
reference. Previously obtained handles (from the old group) keep their prior configuration.
"""
function replace_material(g::CoupledSimulations, member::Symbol, f; domain::Union{Nothing,String} = nothing)
    primaries = getfield(g, :primaries)
    refs = getfield(g, :refs)
    is_primary = haskey(primaries, member)
    is_primary || haskey(refs, member) || throw(ArgumentError("unknown member `$member`"))
    old_sim = is_primary ? primaries[member] : refs[member]
    if domain === nothing
        new_db = replace_material(old_sim.db, f)
    else
        old_sim.db isa DomainBuffers || throw(ArgumentError(
            "`domain` selector requires member `$member` to be a domain dictionary"))
        new_db = replace_material(old_sim.db, domain, f)
    end
    new_sim = Simulation(new_db, old_sim.a, old_sim.aold)
    new_primaries = is_primary ? merge(primaries, NamedTuple{(member,)}((new_sim,))) : primaries
    new_refs = is_primary ? refs : merge(refs, NamedTuple{(member,)}((new_sim,)))
    return CoupledSimulations(new_primaries; refs = new_refs)
end
