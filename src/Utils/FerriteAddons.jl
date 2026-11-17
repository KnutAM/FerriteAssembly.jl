# ================================================= #
# SubDofHandler                                     #
# ================================================= #
_getgrid(sdh::SubDofHandler) = sdh.dh.grid
_getcellset(sdh::SubDofHandler) = sdh.cellset

# Custom functions for SubDofHandler
# Cell coordinates with the type used for the cells in `sdh`
sample_coordinates(sdh::SubDofHandler) = getcoordinates(_getgrid(sdh), first(_getcellset(sdh)))

function create_dofrange(sdh::SubDofHandler)
    return NamedTuple((n => dof_range(sdh, n) for n in Ferrite.getfieldnames(sdh)))
end

# ================================================= #
# Algebraic variables                               #
# ================================================= #
# Ferrite has no public accessors for the following, access to the internal fields is collected here

# The names of all algebraic variables in the dofhandler
_algebraic_names(dh::DofHandler) = dh.algebraic_names
_algebraic_names(sdh::SubDofHandler) = _algebraic_names(sdh.dh)

# The `AlgebraicVariable` with name `name`
function _algebraic_variable(dh::DofHandler, name::Symbol)
    idx = findfirst(==(name), _algebraic_names(dh))
    idx === nothing && throw(ArgumentError(":$name is not an algebraic variable in the DofHandler"))
    return dh.algebraic_variables[idx]
end

# The component number (index into `algebraic_dofs(dh, name)`) of the variable `var` for the
# index `idx`, which is an `Int` (scalars and `Vec`s) or a tuple of `Int`s (tensors).
function _algebraic_component(var::AlgebraicVariable{V}, idx) where {V}
    index = idx isa Integer ? (Int(idx),) : idx
    space = Ferrite._index_space(V)
    if !(index isa Tuple{Vararg{Integer}} && length(index) == length(space) && all(map((i, n) -> 1 ≤ i ≤ n, index, space)))
        throw(ArgumentError("The index $(repr(idx)) is not valid for $(var)"))
    end
    return findfirst(==(Ferrite._canonical_index(V, index)), var.components)::Int
end

# ================================================= #
# reinit!(::NamedTuple, args...)                    #
# ================================================= #
# Should be replaced by MultiCellValues
Ferrite.reinit!(vals::NamedTuple, args...) = map(v->reinit!(v, args...), values(vals))

# ================================================= #
# Ferrite.AbstractAssembler                         #
# ================================================= #
_is_atomic(a::Ferrite.AbstractAssembler) = Ferrite._is_atomic(a)
