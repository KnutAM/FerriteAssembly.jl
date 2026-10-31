using Ferrite.CollectionsOfViews: ArrayOfVectorViews

# `StateVector` stores the state for each entity (cell number) in a domain, indexed
# directly by that number instead of going through a `Dict`'s hashing/bucket lookup.
# `vals` is dense over `1:n`, where `n` is the total number of cells in the grid;
# `set::Vector{Int}` (sorted) holds the numbers that actually belong to this domain and
# defines the documented `AbstractDict`-like interface (`keys`, `values`, `pairs`,
# iteration, `haskey`, `get`, `length`, `==`) - entries in `vals` outside `set` are never
# read through that interface and may be left uninitialized.
#
# When the per-cell state is a `Vector{T}` (one value per quadrature point - the common
# case for material models), `vals::ArrayOfVectorViews{T,1}` flattens every cell's vector
# into one shared buffer instead of allocating a small `Vector` per cell, and
# `getindex`/`get_state` return a `SubArray` view into it rather than a `Vector`. Any
# other per-cell state (a scalar mutable struct, `Nothing`, or a custom `AbstractVector`
# that is not a plain `Vector`, e.g. an immutable NTuple-backed wrapper) is stored as a
# plain dense `Vector{SV}`, unflattened.
#
# The public `getindex`/`setindex!`/`get` below check that `cellnum` belongs to the domain
# (like a `Dict` would) and throw `KeyError` otherwise. The internal, unchecked
# `_raw_getindex`/`_raw_setindex!` (used on the assembly-hot path, where `cellnum` is
# always valid by construction) skip that check - this is exactly the hashing/bucket-
# lookup overhead being removed, so it is kept off the hot path rather than folded into
# the public accessors.
mutable struct StateVector{SV, VV <: AbstractVector{SV}} <: AbstractDict{Int, SV}
    vals::VV
    set::Vector{Int} # sorted
end
@inline _raw_getindex(s::StateVector, cellnum::Int) = s.vals[cellnum]
@inline _raw_setindex!(s::StateVector, v, cellnum::Int) = (s.vals[cellnum] = v; v)
# `ArrayOfVectorViews` has no `setindex!` at the top level (each entry is a fixed-size
# view): replacing a packed cell's state means copying into the existing view instead.
@inline function _raw_setindex!(s::StateVector{SV, <:ArrayOfVectorViews}, v, cellnum::Int) where {SV}
    dst = s.vals[cellnum]
    length(dst) == length(v) || throw(ArgumentError(
        "cannot replace a packed state of length $(length(dst)) with one of length $(length(v))"))
    copyto!(dst, v)
    return v
end
# `fast_getindex` (see Utils/utils.jl) short-circuits a `Nothing`-valued state entirely;
# disambiguate explicitly since `StateVector{Nothing}` also matches that generic method.
fast_getindex(s::StateVector{Nothing}, cellnum) = nothing
fast_getindex(s::StateVector, cellnum) = _raw_getindex(s, Int(cellnum))

function Base.getindex(s::StateVector, cellnum)
    haskey(s, cellnum) || throw(KeyError(cellnum))
    return _raw_getindex(s, Int(cellnum))
end
function Base.setindex!(s::StateVector, v, cellnum)
    haskey(s, cellnum) || throw(KeyError(cellnum))
    return _raw_setindex!(s, v, Int(cellnum))
end
function Base.:(==)(a::StateVector, b::StateVector)
    a.set == b.set || return false
    return all(c -> _raw_getindex(a, c) == _raw_getindex(b, c), a.set)
end
function Base.iterate(s::StateVector, i::Int = 1)
    i > length(s.set) && return nothing
    c = s.set[i]
    return (c => _raw_getindex(s, c), i + 1)
end
Base.length(s::StateVector) = length(s.set)
Base.keys(s::StateVector) = s.set
Base.values(s::StateVector) = (_raw_getindex(s, c) for c in s.set)
Base.pairs(s::StateVector) = (c => _raw_getindex(s, c) for c in s.set)
Base.haskey(s::StateVector, cellnum) = insorted(cellnum, s.set)
Base.get(s::StateVector, cellnum, default) = haskey(s, cellnum) ? _raw_getindex(s, Int(cellnum)) : default

struct StateVariables{SV, VV}
    old::StateVector{SV, VV} # Rule: Referenced during assembly, not changed (ever)
    new::StateVector{SV, VV} # Rule: Updated during assembly, not referenced (before updated)
end

# Internal implementation for `update_states!`, see the docstring on the
# `AbstractDomainBuffer`/`DomainBuffers` method for the meaning of `mode`.
function update_states!(sv::StateVariables; mode::Symbol = :copy)
    if mode === :copy
        _copy_states!(sv.old, sv.new)
    elseif mode === :flip
        _flip_states!(sv)
    else
        throw(ArgumentError("Unknown mode=$(repr(mode)) for update_states!, use :copy or :flip"))
    end
    return sv
end

function _flip_states!(sv::StateVariables)
    tmp = sv.old.vals
    sv.old.vals = sv.new.vals
    sv.new.vals = tmp
    return sv
end

"""
    copy_state(state)

Return a copy of `state` such that the intended mutation of the returned value does not affect `state`.
Used by [`update_states!`](@ref update_states!(::FerriteAssembly.DomainBuffers))'s default
`mode = :copy` and by [`revert_states!`](@ref revert_states!(::FerriteAssembly.DomainBuffers))
to copy one state into the other, when [`copy_state!`](@ref) is not applicable for that value.

For a cell state that is a *mutable* `AbstractArray` (`ismutable(state) == true`), these functions
apply this check element-wise; for any other cell state (including an immutable `AbstractArray`,
e.g. built from `NTuple`s or `StaticArrays`), it applies to the whole state. In both cases, values
for which `isbits(value) == true` are copied by identity internally, without calling `copy_state`
or `copy_state!`. This function has no default method, so it must be overloaded for the type of
any value (the whole cell state, a whole immutable array, or a mutable array's element) for which
`isbits(value) == false`, unless [`copy_state!`](@ref) is overloaded for that value's type instead.
"""
function copy_state end

"""
    copy_state!(dst, src)

Overwrite `dst` in place so that it matches the value of `src`; the return value is ignored.

An optional, allocation-avoiding alternative to [`copy_state`](@ref) for a value that is itself
immutable (so [`copy_state`](@ref) would otherwise need to allocate a full replacement every
time) but wraps a mutable payload that can instead be updated in place, e.g.
`struct MyState; vals::Vector{Float64}; end`. A value's type needs at most one of
`copy_state!(dst, src)` or `copy_state(src)` overloaded — never both — and if
`copy_state!(dst, src)` is applicable, it takes precedence over `copy_state(src)`. Neither is
needed for `isbits` values.
"""
function copy_state! end

@inline function _copy_state_into(dst, src)
    if isbits(src)
        return src
    elseif applicable(copy_state!, dst, src)
        copy_state!(dst, src)
        return dst
    else
        return copy_state(src) # Call user-implementable function
    end
end

function _copy_states!(dst::StateVector, src::StateVector)
    # A `SubArray` view returned from `ArrayOfVectorViews`-packed storage is always safe to
    # mutate element-wise (its parent, the shared flat buffer, is always a mutable `Vector`
    # we own) - unlike the generic branch below, where a value's own `ismutable` must decide
    # (e.g. a *user-returned* `SubArray` whole-cell state must keep going through
    # `copy_state`/`copy_state!` as a whole value, exactly as before this packing existed).
    packed = src.vals isa ArrayOfVectorViews
    for key in src.set
        src_val = _raw_getindex(src, key)
        if isa(src_val, AbstractArray) && (packed || ismutable(src_val))
            dst_val = _raw_getindex(dst, key)
            axes(dst_val) == axes(src_val) || throw(ArgumentError("Dimension mismatch between old and new cell states"))
            @inbounds for i in eachindex(dst_val, src_val)
                src_i = src_val[i]
                if isbits(src_i)
                    dst_val[i] = src_i
                elseif isassigned(dst_val, i) # otherwise nothing meaningful to mutate via copy_state!
                    dst_val[i] = _copy_state_into(dst_val[i], src_i)
                else
                    dst_val[i] = copy_state(src_i)
                end
            end
        elseif isbits(src_val)
            _raw_setindex!(dst, src_val, key)
        else                                                       # `copy_state!`/`copy_state` should
            _raw_setindex!(dst, _copy_state_into(_raw_getindex(dst, key), src_val), key) # be overloaded, not `_copy_state_into`
        end
    end
end

revert_states!(sv::StateVariables) = _copy_states!(sv.new, sv.old)

# Experimental, basically copy!, but use separate name for clarity
function replace_states!(dst::StateVariables, src::StateVariables)
    dst.old.vals = src.old.vals
    dst.new.vals = src.new.vals
    return dst
end

"""
    create_cell_state(material, cellvalues, x, ae, dofrange)

Defaults to returning `nothing`.

Overload this function to create the state which should be passed into the
`element_routine!`/`element_residual!` for the given `material` and `cellvalues`.
`x` is the cell's coordinates, `ae` the element degree of freedom values, and
`dofrange::NamedTuple` containing the local dof range for each field.
As for the element routines, `ae`, is filled with `NaN` unless the global degree
of freedom vector is given to the [`setup_domainbuffer`](@ref) function.
"""
create_cell_state(material, cv, args...) = [nothing for _ in 1:_getnquadpoints(cv)]

_getnquadpoints(fe_v::Ferrite.AbstractValues) = getnquadpoints(fe_v)
_getnquadpoints(nt::NamedTuple) = getnquadpoints(first(nt))

"""
    _create_cell_state(cell, material, cellvalues, a, ae, dofrange, cellnr)

Internal function to reinit and extract the relevant quantities from the
`cell::CellCache`, reinit cellvalues, update `ae` from `a`, and
pass these into the `create_cell_state` function that the user should specify.
"""
function _create_cell_state(coords, dofs, material, cellvalues, a, ae, dofrange, sdh, cellnr)
    getcoordinates!(coords, _getgrid(sdh), cellnr)
    reinit!(cellvalues, getcells(_getgrid(sdh), cellnr), coords)
    celldofs!(dofs, sdh.dh, cellnr)
    _copydofs!(ae, a, dofs)
    return create_cell_state(material, cellvalues, coords, ae, dofrange)
end

# A domain with an empty cellset has no real cell to sample a state type from; used for
# both an empty cell domain and for facet (and other non-cell) domains, which currently
# have no real per-entity state.
_empty_states(ncells::Int) = StateVariables(
    StateVector(Vector{Nothing}(undef, ncells), Int[]), StateVector(Vector{Nothing}(undef, ncells), Int[]))

"""
    create_states(sdh::SubDofHandler, material, cellvalues, a, cellset, dofrange)

Create the `StateVariables` for the cells in
`cellset` (a sorted `Vector{Int}`), where the user should define the
[`create_cell_state`](@ref) function for their `material` (and corresponding
`cellvalues`). `dofrange::NamedTuple` is passed onto `create_cell_state` and contains the
local dof ranges for each field.
"""
function create_states(sdh::SubDofHandler, material, cellvalues, a, cellset, dofrange)
    grid = _getgrid(sdh)
    ncells = getncells(grid)
    isempty(cellset) && return _empty_states(ncells)
    ae = zeros(ndofs_per_cell(sdh))
    coords = getcoordinates(grid, first(cellset))
    dofs = zeros(Int, ndofs_per_cell(sdh))
    # `old`/`new` must be created independently (not aliased); the resulting per-cell
    # `Vector`s are temporary and packed into the domain's shared storage below.
    old_local = [_create_cell_state(coords, dofs, material, cellvalues, a, ae, dofrange, sdh, cellnr) for cellnr in cellset]
    new_local = [_create_cell_state(coords, dofs, material, cellvalues, a, ae, dofrange, sdh, cellnr) for cellnr in cellset]
    return _pack_states(collect(Int, cellset), ncells, old_local, new_local)
end

# Flatten a per-quadpoint `Vector{T}` state into one shared `ArrayOfVectorViews`, avoiding
# one small-`Vector` allocation per cell (the actual point of using `ArrayOfVectorViews`).
function _pack_states(set::Vector{Int}, ncells::Int, old_local::Vector{SV}, new_local::Vector{SV}) where {T, SV <: Vector{T}}
    lengths = zeros(Int, ncells) # 0 for cells outside the domain -> empty view
    for (i, c) in enumerate(set)
        length(old_local[i]) == length(new_local[i]) || throw(ArgumentError(
            "old and new cell state length mismatch for cell $c: $(length(old_local[i])) != $(length(new_local[i]))"))
        lengths[c] = length(old_local[i])
    end
    indices = Vector{Int}(undef, ncells + 1)
    indices[1] = 1
    for c in 1:ncells
        indices[c + 1] = indices[c] + lengths[c]
    end
    ntotal = indices[end] - 1
    old_data = Vector{T}(undef, ntotal)
    new_data = Vector{T}(undef, ntotal)
    for (i, c) in enumerate(set)
        _unsafe_pack_into!(old_data, indices[c], old_local[i])
        _unsafe_pack_into!(new_data, indices[c], new_local[i])
    end
    lin = LinearIndices((ncells,))
    old_vv = ArrayOfVectorViews(indices, old_data, lin)
    new_vv = ArrayOfVectorViews(indices, new_data, lin) # `indices` never mutated: safe to share
    return StateVariables(StateVector(old_vv, set), StateVector(new_vv, set))
end

# Copy `src`'s (possibly partially unassigned, see `create_cell_state` returning
# `Vector{T}(undef, n)`) elements into `dest` starting at `offset`, without ever reading an
# unassigned slot (leaves the corresponding `dest` slot unassigned too, matching a plain
# `Vector{T}(undef, ...)`'s initial state).
function _unsafe_pack_into!(dest::Vector{T}, offset::Int, src::Vector{T}) where {T}
    for i in eachindex(src)
        if isbitstype(T) || isassigned(src, i)
            @inbounds dest[offset + i - 1] = src[i]
        end
    end
    return dest
end

# Any other per-cell state (a scalar mutable struct, `Nothing`, or a custom
# `AbstractVector` that is not a plain `Vector`): dense array, one slot per cell, no
# flattening (this is a `Dict` -> plain-array backing change only).
function _pack_states(set::Vector{Int}, ncells::Int, old_local::Vector{SV}, new_local::Vector{SV}) where {SV}
    old_vals = Vector{SV}(undef, ncells)
    new_vals = Vector{SV}(undef, ncells)
    for (i, c) in enumerate(set)
        old_vals[c] = old_local[i]
        new_vals[c] = new_local[i]
    end
    return StateVariables(StateVector(old_vals, set), StateVector(new_vals, set))
end
