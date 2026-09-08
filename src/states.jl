# Minimal interface for a vector, storage format will probably be updated later. 
mutable struct StateVector{SV}
    vals::Dict{Int, SV}
end
Base.getindex(s::StateVector, cellnum::Int) = s.vals[cellnum]
Base.setindex!(s::StateVector, v, cellnum::Int) = setindex!(s.vals, v, cellnum)
Base.:(==)(a::StateVector, b::StateVector) = (a.vals == b.vals)

struct StateVariables{SV}
    old::StateVector{SV} # Rule: Referenced during assembly, not changed (ever)
    new::StateVector{SV} # Rule: Updated during assembly, not referenced (before updated)
end
StateVariables(old::Dict, new::Dict) = StateVariables(StateVector(old), StateVector(new))

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
    for key in keys(src.vals)
        src_val = src.vals[key]
        if isa(src_val, AbstractArray) && ismutable(src_val)
            dst_val = dst.vals[key]
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
            dst.vals[key] = src_val
        else                                                 # `copy_state!`/`copy_state` should
            dst.vals[key] = _copy_state_into(dst.vals[key], src_val) # be overloaded, not `_copy_state_into`
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

"""
    create_states(sdh::SubDofHandler, material, cellvalues, a, cellset, dofrange)

Create a `Dict` of states for the cells in `cellset`, where the user should 
define the [`create_cell_state`](@ref) function for their `material` (and corresponding `cellvalues`)
`dofrange::NamedTuple` is passed onto `create_cell_state` and contains the local dof ranges for each field. 
"""
function create_states(sdh::SubDofHandler, material, cellvalues, a, cellset, dofrange)
    ae = zeros(ndofs_per_cell(sdh))
    coords = getcoordinates(_getgrid(sdh), first(cellset))
    dofs = zeros(Int, ndofs_per_cell(sdh))
    return Dict(cellnr => _create_cell_state(coords, dofs, material, cellvalues, a, ae, dofrange, sdh, cellnr) for cellnr in cellset)
end
