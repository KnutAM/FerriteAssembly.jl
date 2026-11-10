"""
    algebraic_couplings(buffers...)

Create the coupling descriptors, `Vector{Union{CellCoupling, FacetCoupling}}`, for the domains with
algebraic variables (see the `algebraic_per_cell` keyword of [`DomainSpec`](@ref)) in `buffers`.
Each element in `buffers` can be a domain buffer, a [`DomainBuffers`](@ref FerriteAssembly.DomainBuffers),
or a [`Simulation`](@ref). The output can be passed to Ferrite's `allocate_matrix`, to allocate the
entries required for the algebraic variables in the sparsity pattern, i.e.
```julia
K = allocate_matrix(dh, ch; algebraic_couplings = FerriteAssembly.algebraic_couplings(cellbuffers, facetbuffers))
```
For each domain, all fields in the domain are coupled to all algebraic variables in the domain,
and all algebraic variables in the domain are coupled to each other. If this gives too many entries
in the sparsity pattern, the coupling descriptors can instead be created manually.
"""
function algebraic_couplings(buffers...)
    couplings = Union{CellCoupling, FacetCoupling}[]
    foreach(b -> add_algebraic_couplings!(couplings, b), buffers)
    return couplings
end

add_algebraic_couplings!(couplings, sim::AbstractSimulation) = add_algebraic_couplings!(couplings, get_domainbuffer(sim))
function add_algebraic_couplings!(couplings, dbs::DomainBuffers)
    foreach(db -> add_algebraic_couplings!(couplings, db), values(dbs))
    return couplings
end
function add_algebraic_couplings!(couplings, db::StdDomainBuffer)
    layout = get_layout(get_base(get_itembuffer(db)))
    set = get_set(db)
    (layout === nothing || isempty(set)) && return couplings
    names = layout.names
    pairs = Tuple{Symbol, Symbol}[]
    for fieldname in Ferrite.getfieldnames(db.sdh), name in names
        push!(pairs, (fieldname, name))
    end
    for (i, name) in enumerate(names), j in i:length(names)
        push!(pairs, (name, names[j]))
    end
    push!(couplings, create_coupling(set, pairs))
    return couplings
end

create_coupling(set::Vector{Int}, pairs) = CellCoupling(set; algebraic_coupling = pairs)
create_coupling(set::Vector{FacetIndex}, pairs) = FacetCoupling(set; algebraic_coupling = pairs)
