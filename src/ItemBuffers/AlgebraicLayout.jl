"""
    compute_fieldwise_ranges(fields, attached, standalone)

Compute the field-wise local dof layout of an item (cell or facet) with algebraic dofs.

* `fields::AbstractVector{Pair{Symbol, Int}}`: The spatial fields, `name => ndofs`, in the order of
  the `SubDofHandler`, where `ndofs` is the field's number of dofs per cell.
* `attached::AbstractDict{Symbol, <:AbstractVector{Pair{Symbol, Int}}}`: For a field name, the
  algebraic variables `name => ndofs` whose dofs follow directly after that field's cell dofs.
* `standalone::AbstractVector{Pair{Symbol, Int}}`: The algebraic variables, `name => ndofs`, whose
  dofs follow after all fields.

The local layout is, for each field in order: the field's cell dofs, then the dofs attached to
that field, followed by all standalone variables. Returns a `NamedTuple` with

* `field_sources`: Each field's range in the dofs given by `Ferrite.celldofs`
* `field_destinations`: The first local index of each field's cell dofs
* `field_ranges`: Each field's local range, including the dofs attached to it
* `variable_names` and `variable_ranges`: The names and local ranges of all algebraic variables,
  the attached variables (in field order) before the standalone variables.
* `ndofs`: The total number of local dofs
"""
function compute_fieldwise_ranges(
        fields::AbstractVector{Pair{Symbol, Int}},
        attached::AbstractDict{Symbol, <:AbstractVector{Pair{Symbol, Int}}},
        standalone::AbstractVector{Pair{Symbol, Int}},
    )
    field_sources = UnitRange{Int}[]
    field_destinations = Int[]
    field_ranges = UnitRange{Int}[]
    variable_names = Symbol[]
    variable_ranges = UnitRange{Int}[]
    source_offset = 0
    offset = 0
    for (fieldname, n) in fields
        push!(field_sources, (source_offset + 1):(source_offset + n))
        push!(field_destinations, offset + 1)
        field_start = offset + 1
        source_offset += n
        offset += n
        for (name, nvar) in get(attached, fieldname, Pair{Symbol, Int}[])
            push!(variable_names, name)
            push!(variable_ranges, (offset + 1):(offset + nvar))
            offset += nvar
        end
        push!(field_ranges, field_start:offset)
    end
    for (name, nvar) in standalone
        push!(variable_names, name)
        push!(variable_ranges, (offset + 1):(offset + nvar))
        offset += nvar
    end
    return (; field_sources, field_destinations, field_ranges, variable_names, variable_ranges, ndofs = offset)
end

"""
    AlgebraicLayout

Describes the local dof layout of an item (cell or facet) in a domain that uses algebraic
variables, see the `algebraic_per_cell` keyword of [`DomainSpec`](@ref). The item buffers of
domains without algebraic variables store `nothing` instead.

Internal type, set up by [`setup_domainbuffer`](@ref).
"""
struct AlgebraicLayout{N, NF, V <: NamedTuple}
    ncelldofs::Int                              # Number of dofs given by `Ferrite.celldofs`
    ndofs::Int                                  # Total number of local dofs
    field_sources::NTuple{NF, UnitRange{Int}}   # Each field's range in `Ferrite.celldofs`
    field_destinations::NTuple{NF, Int}         # First local index of each field's cell dofs
    field_ranges::NTuple{NF, UnitRange{Int}}    # Each field's local range
    names::NTuple{N, Symbol}                    # Algebraic variables, in local order
    ranges::NTuple{N, UnitRange{Int}}           # Local range of each algebraic variable
    values::V                                   # NamedTuple{names} with each variable's AlgebraicValues
end

"""
    setup_layout(sdh::SubDofHandler, names::Tuple{Vararg{Symbol}})

Create the [`AlgebraicLayout`](@ref) for a domain on `sdh` where each item uses the algebraic
variables `names`. Returns `nothing` if `names` is empty.
"""
setup_layout(::SubDofHandler, ::Tuple{}) = nothing
function setup_layout(sdh::SubDofHandler, names::Tuple{Vararg{Symbol}})
    dh = sdh.dh
    Ferrite.isclosed(dh) || throw(ArgumentError("The DofHandler must be closed before setting up a domain with algebraic variables"))
    available = _algebraic_names(dh)
    for name in names
        if name ∉ available
            if name ∈ Ferrite.getfieldnames(dh)
                throw(ArgumentError(":$name is a spatial field, but `algebraic_per_cell` must contain algebraic variables"))
            else
                throw(ArgumentError(":$name is not an algebraic variable in the DofHandler (algebraic variables: $available)"))
            end
        end
    end
    allunique(names) || throw(ArgumentError("Each algebraic variable can only be given once per domain, got $names"))

    fields = [fieldname => length(dof_range(sdh, fieldname)) for fieldname in Ferrite.getfieldnames(sdh)]
    standalone = [name => length(algebraic_dofs(dh, name)) for name in names]
    r = compute_fieldwise_ranges(fields, Dict{Symbol, Vector{Pair{Symbol, Int}}}(), standalone)
    values = NamedTuple{names}(map(name -> AlgebraicValues(_algebraic_variable(dh, name)), names))
    return AlgebraicLayout(
        ndofs_per_cell(sdh), r.ndofs,
        Tuple(r.field_sources), Tuple(r.field_destinations), Tuple(r.field_ranges),
        names, Tuple(r.variable_ranges), values)
end

# Number of local dofs for an item
ndofs_per_item(sdh::SubDofHandler, ::Nothing) = ndofs_per_cell(sdh)
ndofs_per_item(::SubDofHandler, layout::AlgebraicLayout) = layout.ndofs

"""
    create_dofrange(sdh::SubDofHandler, layout)

Create the `NamedTuple` with the local dof range for each field in `sdh`, followed by the
range for each algebraic variable in `layout`.
"""
create_dofrange(sdh::SubDofHandler, ::Nothing) = create_dofrange(sdh)
function create_dofrange(sdh::SubDofHandler, layout::AlgebraicLayout)
    fieldnames = Tuple(Ferrite.getfieldnames(sdh))
    return NamedTuple{(fieldnames..., layout.names...)}((layout.field_ranges..., layout.ranges...))
end

"""
    fill_celldofs!(dofs::Vector{Int}, layout, dh, cellnum::Int)

Update the part of the local dofs, `dofs`, that depends on the cell, i.e. the dofs given by
`Ferrite.celldofs`. The dofs of the algebraic variables are not changed, see
[`update_algebraic_dofs!`](@ref).
"""
@inline fill_celldofs!(dofs::Vector{Int}, ::Nothing, dh, cellnum::Int) = celldofs!(dofs, dh, cellnum)
@inline function fill_celldofs!(dofs::Vector{Int}, layout::AlgebraicLayout, dh, cellnum::Int)
    # The dofs from `Ferrite.celldofs` are located first in the local layout
    return celldofs!(view(dofs, 1:layout.ncelldofs), dh, cellnum)
end

"""
    update_algebraic_dofs!(dofs::Vector{Int}, layout, dh)

Write the global dofs of the algebraic variables in `layout` to their positions in the local
dofs, `dofs`. These don't depend on the cell, and are updated once before working a domain.
Nothing is done for a `layout` that is `nothing` (no algebraic variables).
"""
update_algebraic_dofs!(::Vector{Int}, ::Nothing, dh) = nothing
function update_algebraic_dofs!(dofs::Vector{Int}, layout::AlgebraicLayout, dh)
    foreach(layout.names, layout.ranges) do name, range
        copyto!(dofs, first(range), algebraic_dofs(dh, name), 1, length(range))
    end
    return nothing
end

# Copy the dofs of the algebraic variables from one local dof vector to another
function copy_algebraic_dofs!(dst::Vector{Int}, src::Vector{Int}, layout::AlgebraicLayout)
    foreach(range -> copyto!(dst, first(range), src, first(range), length(range)), layout.ranges)
    return nothing
end

@inline _get_algebraic_values(layout::AlgebraicLayout, name::Symbol) = layout.values[name]
_get_algebraic_values(::Nothing, name::Symbol) = throw(ArgumentError(
    "No algebraic variables in the domain, cannot get the algebraic values for :$name"))

# The names of the spatial fields, given all `names` in the dofrange
_spatial_fieldnames(names::Tuple{Vararg{Symbol}}, ::Nothing) = names
function _spatial_fieldnames(names::Tuple{Vararg{Symbol}}, ::AlgebraicLayout{<:Any, NF}) where {NF}
    return ntuple(i -> names[i], Val(NF))
end
