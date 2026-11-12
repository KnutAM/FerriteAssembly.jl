"""
    AlgebraicLayout

Describes the local dof layout of an item (cell or facet) in a domain that uses algebraic
variables, see the `algebraic_per_cell` keyword of [`DomainSpec`](@ref). The local dofs are the
`ncelldofs` dofs given by `Ferrite.celldofs`, followed by the dofs of the algebraic variables
in the order they are given to `DomainSpec`. The item buffers of domains without algebraic
variables store `nothing` instead.

Internal type, set up by [`setup_layout`](@ref).
"""
struct AlgebraicLayout{R <: NamedTuple, V <: NamedTuple}
    ncelldofs::Int  # Number of dofs given by `Ferrite.celldofs`
    ranges::R       # NamedTuple{names} with the local range of each algebraic variable
    values::V       # NamedTuple{names} with the AlgebraicValues of each algebraic variable
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

    ncelldofs = ndofs_per_cell(sdh)
    lengths = map(name -> length(algebraic_dofs(dh, name)), names)
    lasts = ncelldofs .+ cumsum(lengths)
    ranges = NamedTuple{names}(map((l, n) -> (l - n + 1):l, lasts, lengths))
    values = NamedTuple{names}(map(name -> AlgebraicValues(_algebraic_variable(dh, name)), names))
    return AlgebraicLayout(ncelldofs, ranges, values)
end

# Number of local dofs for an item
ndofs_per_item(sdh::SubDofHandler, ::Nothing) = ndofs_per_cell(sdh)
ndofs_per_item(::SubDofHandler, layout::AlgebraicLayout) = layout.ncelldofs + sum(length, layout.ranges)

"""
    create_dofrange(sdh::SubDofHandler, layout)

Create the `NamedTuple` with the local dof range for each field in `sdh`, followed by the
range for each algebraic variable in `layout`.
"""
create_dofrange(sdh::SubDofHandler, ::Nothing) = create_dofrange(sdh)
create_dofrange(sdh::SubDofHandler, layout::AlgebraicLayout) = merge(create_dofrange(sdh), layout.ranges)

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
    foreach(keys(layout.ranges), values(layout.ranges)) do name, range
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
function _spatial_fieldnames(names::NTuple{N, Symbol}, ::AlgebraicLayout{<:NamedTuple{vnames}}) where {N, vnames}
    return ntuple(i -> names[i], Val(N - length(vnames)))
end
