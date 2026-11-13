"""
    AttachAlgebraic(name::Symbol => field::Symbol; indices = nothing)

**Experimental**. Declaration, given in the `algebraic_variables` keyword of [`DomainSpec`](@ref),
that the degrees of freedom of the algebraic variable `name` are attached to the spatial field
`field`. For each item (cell or facet) in the domain, these dofs are then placed directly after the
dofs of `field`, and are included in `dof_range(buffer, field)`. The values object of `field` must
therefore provide one shape function for each attached dof, after the field's ordinary shape functions,
see [Attached variables](@ref).

By default, all components of the variable are attached, in `Ferrite`'s canonical order (which is
the data order of the corresponding `Tensors.jl` type). Alternatively, the attached components can be
given, in the desired order, as `indices`. Each index is an `Int` for scalar and `Vec` variables,
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
order of the `SubDofHandler`), the local dofs are the field's dofs given by `Ferrite.celldofs`, followed
by the dofs of the variables attached to the field (see [`AttachAlgebraic`](@ref FerriteAssembly.AttachAlgebraic)),
in declaration order. These are followed by the dofs of the standalone variables, in declaration order.
The item buffers of domains without algebraic variables store `nothing` instead.

Internal type, set up by [`setup_layout`](@ref).
"""
struct AlgebraicLayout{R <: NamedTuple, C <: NamedTuple, V <: NamedTuple, FR <: NamedTuple, OR <: NamedTuple, F}
    ncelldofs::Int      # Number of dofs given by `Ferrite.celldofs`
    ranges::R           # NamedTuple{names} with the local range of each algebraic variable
    components::C       # NamedTuple{names} with the components (indices into `algebraic_dofs(dh, name)`) of each variable
    values::V           # NamedTuple{names} with the AlgebraicValues of each variable (`nothing` if not all components are used)
    field_ranges::FR    # NamedTuple{fieldnames} with the local range of each field, including attached dofs
    own_ranges::OR      # NamedTuple with the local range of each field's own dofs, followed by `ranges`
    fields::F           # `nothing` if no variables are attached, otherwise a tuple with the
                        # range in `celldofs` and the first local dof for each field's own dofs
end

"""
    setup_layout(sdh::SubDofHandler, declarations::Tuple)

Create the [`AlgebraicLayout`](@ref) for a domain on `sdh` with the algebraic variables
in `declarations`, given as the `algebraic_variables` keyword to [`DomainSpec`](@ref).
Returns `nothing` if `declarations` is empty.
"""
setup_layout(::SubDofHandler, ::Tuple{}) = nothing
function setup_layout(sdh::SubDofHandler, declarations::Tuple)
    dh = sdh.dh
    Ferrite.isclosed(dh) || throw(ArgumentError("The DofHandler must be closed before setting up a domain with algebraic variables"))
    names = map(_declared_name, declarations)
    available = _algebraic_names(dh)
    fieldnames = Tuple(Ferrite.getfieldnames(sdh))
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
    components = NamedTuple{names}(map(decl -> _attached_components(dh, decl), declarations))

    # Local ranges, see the docstring of `AlgebraicLayout`
    last_dof = 0
    var_ranges = Dict{Symbol, UnitRange{Int}}()
    field_ranges = Dict{Symbol, UnitRange{Int}}()
    own_field_ranges = Dict{Symbol, UnitRange{Int}}()
    for field in fieldnames
        own_field_ranges[field] = (last_dof + 1):(last_dof + length(dof_range(sdh, field)))
        last_dof = last(own_field_ranges[field])
        for decl in declarations
            if decl isa AttachAlgebraic && decl.field === field
                var_ranges[decl.name] = (last_dof + 1):(last_dof + length(components[decl.name]))
                last_dof = last(var_ranges[decl.name])
            end
        end
        field_ranges[field] = first(own_field_ranges[field]):last_dof
    end
    for decl in declarations
        if decl isa Symbol
            var_ranges[decl] = (last_dof + 1):(last_dof + length(components[decl]))
            last_dof = last(var_ranges[decl])
        end
    end
    ranges = NamedTuple{names}(map(name -> var_ranges[name], names))
    values = NamedTuple{names}(map(name -> _layout_values(dh, name, components[name]), names))
    any_attached = any(decl -> decl isa AttachAlgebraic, declarations)
    fields = any_attached ? map(f -> (dof_range(sdh, f), first(own_field_ranges[f])), fieldnames) : nothing
    return AlgebraicLayout(
        ndofs_per_cell(sdh), ranges, components, values,
        NamedTuple{fieldnames}(map(f -> field_ranges[f], fieldnames)),
        merge(NamedTuple{fieldnames}(map(f -> own_field_ranges[f], fieldnames)), ranges),
        fields)
end

# The components of the variable used in the domain, as indices into `algebraic_dofs(dh, name)`
function _attached_components(dh::DofHandler, name::Symbol)
    return ntuple(identity, Ferrite.n_algebraic_dofs(_algebraic_variable(dh, name)))
end
_attached_components(dh::DofHandler, decl::AttachAlgebraic{Nothing}) = _attached_components(dh, decl.name)
function _attached_components(dh::DofHandler, decl::AttachAlgebraic{<:Tuple})
    isempty(decl.indices) && throw(ArgumentError("The `indices` for :$(decl.name) cannot be empty"))
    var = _algebraic_variable(dh, decl.name)
    components = map(idx -> _algebraic_component(var, idx), decl.indices)
    if !allunique(components)
        throw(ArgumentError("The `indices`, $(decl.indices), for :$(decl.name) refer to the same component more than once"))
    end
    return components
end

# AlgebraicValues if the whole variable is used, otherwise `nothing`
function _layout_values(dh::DofHandler, name::Symbol, components::Tuple)
    var = _algebraic_variable(dh, name)
    return components == ntuple(identity, Ferrite.n_algebraic_dofs(var)) ? AlgebraicValues(var) : nothing
end

# Number of local dofs for an item
ndofs_per_item(sdh::SubDofHandler, ::Nothing) = ndofs_per_cell(sdh)
ndofs_per_item(::SubDofHandler, layout::AlgebraicLayout) = layout.ncelldofs + sum(length, layout.ranges)

# Length of the scratch vector required by `fill_celldofs!`
scratch_length(::Nothing) = 0
scratch_length(layout::AlgebraicLayout) = layout.fields === nothing ? 0 : layout.ncelldofs

"""
    create_dofrange(sdh::SubDofHandler, layout)

Create the `NamedTuple` with the local dof range for each field in `sdh` (including attached
dofs), followed by the range for each algebraic variable in `layout`.
"""
create_dofrange(sdh::SubDofHandler, ::Nothing) = create_dofrange(sdh)
create_dofrange(::SubDofHandler, layout::AlgebraicLayout) = merge(layout.field_ranges, layout.ranges)

"""
    fill_celldofs!(dofs::Vector{Int}, scratch::Vector{Int}, layout, dh, cellnum::Int)

Update the part of the local dofs, `dofs`, that depends on the cell, i.e. the dofs given by
`Ferrite.celldofs`. The dofs of the algebraic variables are not changed, see
[`update_algebraic_dofs!`](@ref). If variables are attached to fields, `scratch` is used to
get the cell dofs before copying each field's dofs to its local position.
"""
@inline fill_celldofs!(dofs::Vector{Int}, ::Vector{Int}, ::Nothing, dh, cellnum::Int) = celldofs!(dofs, dh, cellnum)
@inline function fill_celldofs!(dofs::Vector{Int}, scratch::Vector{Int}, layout::AlgebraicLayout, dh, cellnum::Int)
    return _fill_celldofs!(dofs, scratch, layout.fields, layout.ncelldofs, dh, cellnum)
end
@inline function _fill_celldofs!(dofs::Vector{Int}, ::Vector{Int}, ::Nothing, ncelldofs::Int, dh, cellnum::Int)
    # Without attached variables, the dofs from `Ferrite.celldofs` are located first in the local layout
    return celldofs!(view(dofs, 1:ncelldofs), dh, cellnum)
end
@inline function _fill_celldofs!(dofs::Vector{Int}, scratch::Vector{Int}, fields::Tuple, ::Int, dh, cellnum::Int)
    celldofs!(scratch, dh, cellnum)
    for (src, dst) in fields
        copyto!(dofs, dst, scratch, first(src), length(src))
    end
    return dofs
end

"""
    update_algebraic_dofs!(dofs::Vector{Int}, layout, dh)

Write the global dofs of the algebraic variables in `layout` to their positions in the local
dofs, `dofs`. These don't depend on the cell, and are updated once before working a domain.
Nothing is done for a `layout` that is `nothing` (no algebraic variables).
"""
update_algebraic_dofs!(::Vector{Int}, ::Nothing, dh) = nothing
function update_algebraic_dofs!(dofs::Vector{Int}, layout::AlgebraicLayout, dh)
    foreach(keys(layout.ranges), values(layout.ranges), values(layout.components)) do name, range, components
        global_dofs = algebraic_dofs(dh, name)
        for (i, c) in enumerate(components)
            dofs[first(range) + i - 1] = global_dofs[c]
        end
    end
    return nothing
end

# Copy the dofs of the algebraic variables from one local dof vector to another
function copy_algebraic_dofs!(dst::Vector{Int}, src::Vector{Int}, layout::AlgebraicLayout)
    foreach(range -> copyto!(dst, first(range), src, first(range), length(range)), layout.ranges)
    return nothing
end

@inline _get_algebraic_values(layout::AlgebraicLayout, name::Symbol) = _check_algebraic_values(layout.values[name], name)
_get_algebraic_values(::Nothing, name::Symbol) = throw(ArgumentError(
    "No algebraic variables in the domain, cannot get the algebraic values for :$name"))
@inline _check_algebraic_values(av::AlgebraicValues, ::Symbol) = av
_check_algebraic_values(::Nothing, name::Symbol) = throw(ArgumentError(
    "Only some components of :$name are used in the domain, cannot get the algebraic values"))

@inline _own_ranges(::Nothing, dofrange::NamedTuple) = dofrange
@inline _own_ranges(layout::AlgebraicLayout, ::NamedTuple) = layout.own_ranges

# The names of the spatial fields, given all `names` in the dofrange
_spatial_fieldnames(names::Tuple{Vararg{Symbol}}, ::Nothing) = names
function _spatial_fieldnames(names::NTuple{N, Symbol}, ::AlgebraicLayout{<:NamedTuple{vnames}}) where {N, vnames}
    return ntuple(i -> names[i], Val(N - length(vnames)))
end
