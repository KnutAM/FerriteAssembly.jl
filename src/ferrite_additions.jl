# This files contains function that could be PR:ed into Ferrite.jl,
# but hasn't yet (either due to time or because it makes sense to 
# test the usability a bit first)

@inline _getgrid(dh::Union{DofHandler,MixedDofHandler}) = dh.grid

Ferrite.ndofs_per_cell(dh::MixedDofHandler, fh::FieldHandler) = Ferrite.ndofs_per_cell(dh, first(fh.cellset))
Ferrite.nnodes_per_cell(dh::MixedDofHandler, fh::FieldHandler) = Ferrite.nnodes_per_cell(dh, first(fh.cellset))
Ferrite.nnodes_per_cell(dh::DofHandler) = Ferrite.nnodes_per_cell(_getgrid(dh), 1)

Ferrite.getdim(::Union{MixedDofHandler{dim}, DofHandler{dim}}) where dim = dim
Ferrite.getncells(dh::Union{DofHandler,MixedDofHandler}) = getncells(_getgrid(dh))

Ferrite.getcellset(dh::Union{DofHandler,MixedDofHandler}, args...) = getcellset(_getgrid(dh), args...)

# If a tuple or named tuple is given to reinit! - try to reinit each argument
Ferrite.reinit!(vals::NTuple{N,Ferrite.Values}, args...) where N = foreach(v->reinit!(v, args...), vals) 
Ferrite.reinit!(vals::NamedTuple, args...) = foreach(v->reinit!(v, args...), vals)


# Warning: Black magic up ahead...
# The above solution is very broad in its overload, 
# and representing the cellvalues for multiple fields 
# as a tuple or named tuple can become disambigeous in 
# some cases, see e.g. creation of CellBuffer. 
# An alternative solution, currently not used, but left 
# here as notes, is to introduce a `MultiValues` type.
# It works just like a named tuple, e.g. mv.u and mv[:u] 
# works to access the field :u.
# But it has a well-defined type that requires 
# 1) All elements must be <:Ferrite.Values
# 2) All elements must have the same supertype, i.e. 
#    - Be all Cell*Values or Face*Values 
#    - Have same `dim`, `T`, and `RefShape` type parameters. 
# 
# And it has a nice constructor, 
# `MultiValues(p=CellScalarValues(...), u=CellVectorValues(...))`
#
# Note: It is not possible to access the `vals` field of `mv::MultiValues`
# as mv.vals, this must be done by `getfield(mv, :vals)` as the dot 
# notation refers to the fields of the `vals` tuple!

#=
struct MultiValues{V<:NamedTuple}
    vals::V
    function MultiValues(vals::V) where V
        VT = supertype(typeof(first(vals)))
        isa(first(vals), Ferrite.Values) || throw(ArgumentError("values must be <:Ferrite.Values"))
        all(v->isa(v,VT), vals) || throw(ArgumentError("Incompatible `Values`: $vals"))
        return new{T}(vals)
    end
end
MultiValues(;kwargs...) = MultiValues(NamedTuple(kwargs))

Base.getproperty(mv::MultiValues, f::Symbol) = getproperty(getfield(mv,:vals), f)
Base.getindex(mv::MultiValues, f::Symbol) = getindex(getfield(mv,:vals), f)
Ferrite.reinit!(mv::MultiValues, args...) = foreach(v->reinit!(v, args...), getfield(mv,:vals))

=#