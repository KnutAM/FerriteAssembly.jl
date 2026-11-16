"""
    algebraic_couplings(buffers...; exclude_attached::Bool = false)

Create a `Vector` of coupling descriptors for the domains with algebraic variables (see the
`algebraic_variables` keyword of [`DomainSpec`](@ref)) in `buffers`.
Each element in `buffers` can be a domain buffer, a [`DomainBuffers`](@ref FerriteAssembly.DomainBuffers),
or a [`Simulation`](@ref). The output can be passed to Ferrite's `allocate_matrix`, to allocate the
entries required for the algebraic variables in the sparsity pattern, i.e.
```julia
algebraic_couplings = FerriteAssembly.algebraic_couplings(cellbuffers, facetbuffers)
K = allocate_matrix(dh, ch; algebraic_couplings)
```
For each domain, all fields in the domain are coupled to all algebraic variables in the domain,
and all algebraic variables in the domain are coupled to each other. This includes all components of
a variable, also those that are not attached in the domain, see [`AttachAlgebraic`](@ref FerriteAssembly.AttachAlgebraic).
If more coupling entries are desired, they can be `push!`ed to `algebraic_couplings`.

If `exclude_attached = true`, the variables that are attached to a field are excluded, such that only the
standalone variables are coupled. The entries for the attached variables must then be added with
[`add_algebraic_entries!`](@ref FerriteAssembly.add_algebraic_entries!), which only adds the entries for the
attached components.
"""
function algebraic_couplings(buffers...; exclude_attached::Bool = false)
    couplings = Union{CellCoupling, FacetCoupling, AlgebraicCoupling}[]
    foreach(b -> add_algebraic_couplings!(couplings, b, exclude_attached), buffers)
    return couplings
end

add_algebraic_couplings!(couplings, sim::AbstractSimulation, exclude_attached::Bool) = add_algebraic_couplings!(couplings, get_domainbuffer(sim), exclude_attached)
function add_algebraic_couplings!(couplings, dbs::DomainBuffers, exclude_attached::Bool)
    foreach(db -> add_algebraic_couplings!(couplings, db, exclude_attached), values(dbs))
    return couplings
end
function add_algebraic_couplings!(couplings, db::StdDomainBuffer, exclude_attached::Bool)
    layout = get_layout(get_base(get_itembuffer(db)))
    set = get_set(db)
    (layout === nothing || isempty(set)) && return couplings
    names = exclude_attached ? keys(layout.standalone_values) : _variable_names(layout)
    isempty(names) && return couplings
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

"""
    FerriteAssembly.add_algebraic_entries!(sp::Ferrite.AbstractSparsityPattern, buffers...)

**Experimental**. Add the sparsity entries for the attached components of the algebraic variables
(see [`AttachAlgebraic`](@ref FerriteAssembly.AttachAlgebraic)) in the domains in `buffers` to `sp`.
As for [`algebraic_couplings`](@ref FerriteAssembly.algebraic_couplings), each element in `buffers` can be a
domain buffer, a [`DomainBuffers`](@ref FerriteAssembly.DomainBuffers), or a [`Simulation`](@ref).
For each item in a domain, the attached components are coupled to all degrees of freedom of the (adjacent) cell,
to each other, and to the standalone algebraic variables in the domain. In contrast to `algebraic_couplings`,
the components of a variable that are not attached in any domain don't get entries, instead of dense rows and
columns.

The entries for the fields, the standalone variables, and the constraints are not added. Hence, the
full pattern is created by
```julia
sp = init_sparsity_pattern(dh)
algebraic_couplings = FerriteAssembly.algebraic_couplings(buffers; exclude_attached = true)
add_sparsity_entries!(sp, dh; algebraic_couplings)  # Entries for the fields and standalone variables
FerriteAssembly.add_algebraic_entries!(sp, buffers)
add_constraint_entries!(sp, ch)                     # Must be added last
K = allocate_matrix(sp)
```
"""
function add_algebraic_entries!(sp::Ferrite.AbstractSparsityPattern, buffers...)
    # The rows of the attached dofs can become dense. Their columns are therefore collected
    # for all items, and added in sorted order at the end.
    collected_cols = Dict{Int, Vector{Int}}()
    foreach(b -> _add_algebraic_entries!(sp, collected_cols, b), buffers)
    for (row, cols) in collected_cols
        for col in unique!(sort!(cols))
            Ferrite.add_entry!(sp, row, col)
        end
    end
    return sp
end

_add_algebraic_entries!(sp, collected_cols, sim::AbstractSimulation) = _add_algebraic_entries!(sp, collected_cols, get_domainbuffer(sim))
function _add_algebraic_entries!(sp, collected_cols, dbs::DomainBuffers)
    foreach(db -> _add_algebraic_entries!(sp, collected_cols, db), values(dbs))
    return nothing
end
function _add_algebraic_entries!(sp, collected_cols, db::StdDomainBuffer)
    layout = get_layout(get_base(get_itembuffer(db)))
    set = get_set(db)
    (layout === nothing || isempty(set)) && return nothing
    dh = db.sdh.dh
    # Global dofs of the attached components (attached_dofs) and of the standalone variables (standalone_dofs)
    attached_dofs = Int[]
    foreach(keys(layout.attached), values(layout.attached)) do name, (; components)
        global_dofs = algebraic_dofs(dh, name)
        foreach(c -> push!(attached_dofs, global_dofs[c]), components)
    end
    isempty(attached_dofs) && return nothing
    standalone_dofs = mapreduce(name -> algebraic_dofs(dh, name), vcat, keys(layout.standalone_values); init = Int[])
    for row in attached_dofs
        append!(get!(Vector{Int}, collected_cols, row), attached_dofs, standalone_dofs)
    end
    for row in standalone_dofs, col in attached_dofs
        Ferrite.add_entry!(sp, row, col)
    end
    cdofs = Int[]
    for cellnr in _item_cells(set)
        resize!(cdofs, ndofs_per_cell(dh, cellnr))
        celldofs!(cdofs, dh, cellnr)
        for row in attached_dofs
            append!(collected_cols[row], cdofs)
        end
        for row in cdofs, col in attached_dofs
            Ferrite.add_entry!(sp, row, col)
        end
    end
    return nothing
end

# The (unique) cells of the items in a domain's set
_item_cells(set::Vector{Int}) = set
_item_cells(set::Vector{FacetIndex}) = unique!([facet[1] for facet in set])
