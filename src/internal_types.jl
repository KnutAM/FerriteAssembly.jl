"""
    FixedDict{K,V} <: AbstractDict{K,V}

A `Dict`-backed `AbstractDict` that is publicly fixed (immutable) but internally mutable when
safe. `getindex`, `haskey`, `get`, `keys`, `values`, `pairs`, `length`, `isempty` and
iteration forward to the wrapped `Dict{K,V}`, so their performance and allocation behavior
are identical to a plain `Dict`; no `setindex!` is defined, and `getproperty` always throws,
so nothing outside this file can reach or mutate the wrapped `Dict`.

Internally, the wrapped `Dict` remains reachable via `getfield(d, :dict)`, and can still be
mutated in place when that is actually safe to do - see `_internal_setindex!`.
"""
struct FixedDict{K,V} <: AbstractDict{K,V}
    dict::Dict{K,V}
    global _unsafe_fixeddict(d::Dict{K,V}) where {K,V} = new{K,V}(d)
end
FixedDict(itr) = _unsafe_fixeddict(Dict(itr))
FixedDict{K,V}() where {K,V} = _unsafe_fixeddict(Dict{K,V}())

Base.getproperty(::FixedDict, ::Symbol) = error(
    "`FixedDict` fields are not accessible via `.`; this type has no public fields")
Base.propertynames(::FixedDict) = ()

Base.length(d::FixedDict) = length(getfield(d, :dict))
Base.isempty(d::FixedDict) = isempty(getfield(d, :dict))
Base.haskey(d::FixedDict, k) = haskey(getfield(d, :dict), k)
Base.get(d::FixedDict, k, default) = get(getfield(d, :dict), k, default)
Base.getindex(d::FixedDict, k) = getfield(d, :dict)[k]
Base.iterate(d::FixedDict) = iterate(getfield(d, :dict))
Base.iterate(d::FixedDict, state) = iterate(getfield(d, :dict), state)
# `keys`/`values`/`pairs` are intentionally NOT forwarded to the wrapped `Dict`: doing so
# (e.g. `Base.keys(d::FixedDict) = keys(getfield(d, :dict))`) would return a `KeySet`/
# `ValueIterator` wrapping the *raw* `Dict`, whose own `.dict` field (unprotected) would hand
# it right back out. Falling through to `AbstractDict`'s generic defaults instead wraps `d`
# itself (so `keys(d).dict`/`values(d).dict` just returns `d`, which is no less protected) -
# at no performance or allocation cost, since those defaults only need `length`/`iterate`.

"""
    _internal_setindex!(d::FixedDict, kv::Pair)

Insert `kv` into `d`'s wrapped `Dict` in place (replacing, not shadowing, an already-present
key - this is a real `Dict`, not `Base.ImmutableDict`), bypassing the public immutability.

**Only call this on a `FixedDict` that is not part of its owning type's own *documented*
public interface.** `Simulation.db` is such a documented interface (`Simulation`'s own
docstring advertises `.db` as something user code reads directly) and must never be mutated
this way - any code holding a `Simulation` could be relying on its `.db` staying exactly as it
was. `LoadHandler.nbcs`/`.bodyloads`, by contrast, are undocumented implementation-detail
fields: `LoadHandler`'s own documented interface is `add!`/`apply!` only, so mutating them in
place here is the same category of "deliberately bypassing the documented interface" as
reading `getfield(d, :dict)` on a `FixedDict` itself - not something this package
protects against. When calling this on a type whose `FixedDict` field *is* public and
documented, build a full replacement with [`FixedDict`](@ref) instead and reassign the
field, the way `replace_material` does.
"""
_internal_setindex!(d::FixedDict, kv::Pair) = (getfield(d, :dict)[kv.first] = kv.second; d)
