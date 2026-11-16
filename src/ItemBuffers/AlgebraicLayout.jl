"""
    AttachAlgebraic(name::Symbol => field::Symbol; indices = nothing)

**Experimental**. Declaration, given in the `algebraic_variables` keyword of [`DomainSpec`](@ref),
that the degrees of freedom of the algebraic variable `name` are attached to the spatial field
`field`. For each item (cell or facet) in the domain, these dofs are then placed directly after the
dofs of `field`, and are included in `dof_range(buffer, field)`. The values object of `field` must
therefore provide one shape function for each attached dof, after the field's ordinary shape functions,
see [Attached variables](@ref). The attached variable is only seen through `field`, i.e. `name` is not
available in the item buffer (e.g. for `dof_range` and `get_algebraic_values`), nor in the `dofrange`
given to [`create_cell_state`](@ref).

The attached components are given as `indices`, in the order of the attached shape functions.
`indices` can only be omitted if the variable has a single component, e.g. a scalar. Each index is an `Int` for scalar and `Vec` variables,
and a tuple, e.g. `(3, 3)`, for tensorial variables. For symmetric tensors, `(i, j)` and `(j, i)`
refer to the same component. Components that are not attached in any domain don't get any contributions,
and must be constrained.
"""
struct AttachAlgebraic{I}
    name::Symbol
    field::Symbol
    indices::I # `nothing` (all components) or a tuple of indices
end
function AttachAlgebraic(p::Pair{Symbol, Symbol}; indices = nothing)
    if !(indices === nothing || indices isa Tuple)
        throw(ArgumentError("`indices` must be a tuple of indices or `nothing`, got $(repr(indices))"))
    end
    return AttachAlgebraic(p.first, p.second, indices)
end

_declared_name(name::Symbol) = name
_declared_name(a::AttachAlgebraic) = a.name

"""
    AlgebraicLayout

Describes the local dof layout of an item (cell or facet) in a domain that uses algebraic
variables, see the `algebraic_variables` keyword of [`DomainSpec`](@ref). For each field (in the
order of their dofs in `Ferrite.celldofs`), the local dofs are the field's dofs given by `Ferrite.celldofs`, followed
by the dofs of the variables attached to the field (see [`AttachAlgebraic`](@ref FerriteAssembly.AttachAlgebraic)),
in declaration order. These are followed by the dofs of the standalone variables, in declaration order.
Only the fields (with their attached dofs) and the standalone variables are visible in the buffer's `dofrange`,
which is stored in the [`DofBuffer`](@ref FerriteAssembly.DofBuffer).
The item buffers of domains without algebraic variables store `nothing` instead.

Internal type, set up by [`setup_layout`](@ref).
"""
struct AlgebraicLayout{OR <: NamedTuple, AC <: NamedTuple, SV <: NamedTuple}
    # Relation to `Ferrite.celldofs`
    ncelldofs::Int              # Number of dofs given by `Ferrite.celldofs`
    own_ranges::OR              # Local range of each field's own dofs (those given by `celldofs`), followed by the range
                                # of each algebraic variable. Used for the residual scaling and to fill the cell dofs.
    # Attached variables
    attached_components::AC     # name => the attached components (indices into `algebraic_dofs(dh, name)`)
    # Standalone variables
    standalone_values::SV       # name => AlgebraicValues
end

"""
    setup_layout(sdh::SubDofHandler, declarations::Tuple) -> (layout, dofrange)

Create the [`AlgebraicLayout`](@ref) for a domain on `sdh` with the algebraic variables
in `declarations`, given as the `algebraic_variables` keyword to [`DomainSpec`](@ref), and the
`NamedTuple` with the local dof range for each field (including attached dofs), followed by the
range for each standalone variable. If `declarations` is empty, `layout` is `nothing`.
"""
setup_layout(sdh::SubDofHandler, ::Tuple{}) = (nothing, create_dofrange(sdh))
function setup_layout(sdh::SubDofHandler, declarations::Tuple)
    dh = sdh.dh
    Ferrite.isclosed(dh) || throw(ArgumentError("The DofHandler must be closed before setting up a domain with algebraic variables"))
    names = map(_declared_name, declarations)
    available = _algebraic_names(dh)
    # Sorted by the position in `celldofs`, such that the fields' own dofs are in the same order locally
    fieldnames = Tuple(sort!(collect(Ferrite.getfieldnames(sdh)); by = f -> first(dof_range(sdh, f))))
    for name in names
        if name ∉ available
            if name ∈ Ferrite.getfieldnames(dh)
                throw(ArgumentError(":$name is a spatial field, but `algebraic_variables` must contain algebraic variables"))
            else
                throw(ArgumentError(":$name is not an algebraic variable in the DofHandler (algebraic variables: $available)"))
            end
        end
    end
    allunique(names) || throw(ArgumentError("Each algebraic variable can only be given once per domain, got $names"))
    for decl in declarations
        if decl isa AttachAlgebraic && decl.field ∉ fieldnames
            throw(ArgumentError("Cannot attach :$(decl.name) to :$(decl.field), which is not a field in the SubDofHandler (fields: $fieldnames)"))
        end
    end
    attached = Tuple(decl for decl in declarations if decl isa AttachAlgebraic)
    standalone_names = Tuple(decl for decl in declarations if decl isa Symbol)
    attached_components = NamedTuple{map(_declared_name, attached)}(map(decl -> _attached_components(dh, decl), attached))

    # Local ranges, see the docstring of `AlgebraicLayout`
    last_dof = 0
    own_ranges = Dict{Symbol, UnitRange{Int}}()
    field_ranges = Dict{Symbol, UnitRange{Int}}()
    for field in fieldnames
        own_ranges[field] = (last_dof + 1):(last_dof + length(dof_range(sdh, field)))
        last_dof = last(own_ranges[field])
        for decl in attached
            if decl.field === field
                own_ranges[decl.name] = (last_dof + 1):(last_dof + length(attached_components[decl.name]))
                last_dof = last(own_ranges[decl.name])
            end
        end
        field_ranges[field] = first(own_ranges[field]):last_dof
    end
    for name in standalone_names
        own_ranges[name] = (last_dof + 1):(last_dof + length(algebraic_dofs(dh, name)))
        last_dof = last(own_ranges[name])
    end
    standalone_values = NamedTuple{standalone_names}(map(name -> AlgebraicValues(_algebraic_variable(dh, name)), standalone_names))
    layout = AlgebraicLayout(
        ndofs_per_cell(sdh), NamedTuple{(fieldnames..., names...)}(map(name -> own_ranges[name], (fieldnames..., names...))),
        attached_components, standalone_values)
    dofrange = merge(NamedTuple{fieldnames}(map(f -> field_ranges[f], fieldnames)),
        NamedTuple{standalone_names}(map(name -> own_ranges[name], standalone_names)))
    return layout, dofrange
end

# The attached components of the variable, as indices into `algebraic_dofs(dh, name)`
function _attached_components(dh::DofHandler, decl::AttachAlgebraic{Nothing})
    n = length(algebraic_dofs(dh, decl.name))
    if n != 1
        throw(ArgumentError(":$(decl.name) has $n components, so the attached components must be given as `indices`, in the order of the attached shape functions"))
    end
    return (1,)
end
function _attached_components(dh::DofHandler, decl::AttachAlgebraic{<:Tuple})
    isempty(decl.indices) && throw(ArgumentError("The `indices` for :$(decl.name) cannot be empty"))
    var = _algebraic_variable(dh, decl.name)
    components = map(idx -> _algebraic_component(var, idx), decl.indices)
    if !allunique(components)
        throw(ArgumentError("The `indices`, $(decl.indices), for :$(decl.name) refer to the same component more than once"))
    end
    return components
end

# Number of local dofs for an item
ndofs_per_item(sdh::SubDofHandler, ::Nothing) = ndofs_per_cell(sdh)
ndofs_per_item(::SubDofHandler, layout::AlgebraicLayout) = sum(length, layout.own_ranges)

# The names of the fields and of the algebraic variables in the domain
function _fieldnames(::AlgebraicLayout{<:NamedTuple{K}, <:NamedTuple{A}, <:NamedTuple{S}}) where {K, A, S}
    return ntuple(i -> K[i], Val(length(K) - length(A) - length(S)))
end
_variable_names(layout::AlgebraicLayout) = (keys(layout.attached_components)..., keys(layout.standalone_values)...)

# The part of the local `dofs` that is given by `Ferrite.celldofs`: `dofs` itself without algebraic variables,
# otherwise a view of the fields' own dofs
field_dofs_view(dofs::Vector{Int}, ::Nothing) = dofs
function field_dofs_view(dofs::Vector{Int}, layout::AlgebraicLayout)
    indices = mapreduce(f -> collect(layout.own_ranges[f]), vcat, _fieldnames(layout); init = Int[])
    return indices == 1:layout.ncelldofs ? view(dofs, 1:layout.ncelldofs) : view(dofs, indices)
end

"""
    update_algebraic_dofs!(dofs::Vector{Int}, layout, dh)

Write the global dofs of the algebraic variables in `layout` to their positions in the local
dofs, `dofs`. These don't depend on the cell, and are updated once before working a domain.
Nothing is done for a `layout` that is `nothing` (no algebraic variables).
"""
update_algebraic_dofs!(::Vector{Int}, ::Nothing, dh) = nothing
function update_algebraic_dofs!(dofs::Vector{Int}, layout::AlgebraicLayout, dh)
    foreach(keys(layout.attached_components), values(layout.attached_components)) do name, components
        global_dofs = algebraic_dofs(dh, name)
        for (i, c) in enumerate(components)
            dofs[first(layout.own_ranges[name]) + i - 1] = global_dofs[c]
        end
    end
    foreach(keys(layout.standalone_values)) do name
        global_dofs = algebraic_dofs(dh, name)
        copyto!(dofs, first(layout.own_ranges[name]), global_dofs, 1, length(global_dofs))
    end
    return nothing
end

# Copy the dofs of the algebraic variables from one local dof vector to another
function copy_algebraic_dofs!(dst::Vector{Int}, src::Vector{Int}, layout::AlgebraicLayout)
    foreach(_variable_names(layout)) do name
        range = layout.own_ranges[name]
        copyto!(dst, first(range), src, first(range), length(range))
    end
    return nothing
end

@inline function _get_algebraic_values(layout::AlgebraicLayout, name::Symbol)
    haskey(layout.standalone_values, name) || _throw_not_standalone(layout, name)
    return layout.standalone_values[name]
end
_get_algebraic_values(::Nothing, name::Symbol) = throw(ArgumentError(
    "No algebraic variables in the domain, cannot get the algebraic values for :$name"))
@noinline function _throw_not_standalone(layout::AlgebraicLayout, name::Symbol)
    if haskey(layout.attached_components, name)
        throw(ArgumentError(":$name is attached to a field, and can only be accessed via that field"))
    end
    throw(ArgumentError(":$name is not a standalone algebraic variable in the domain (standalone variables: $(keys(layout.standalone_values)))"))
end

@inline _scaling_ranges(::Nothing, dofrange::NamedTuple) = dofrange
@inline _scaling_ranges(layout::AlgebraicLayout, ::NamedTuple) = layout.own_ranges

# The names of the spatial fields, given all `names` in the dofrange
_spatial_fieldnames(names::Tuple{Vararg{Symbol}}, ::Nothing) = names
_spatial_fieldnames(::Tuple{Vararg{Symbol}}, layout::AlgebraicLayout) = _fieldnames(layout)
