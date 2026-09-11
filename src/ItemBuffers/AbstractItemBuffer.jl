abstract type AbstractItemBuffer end

# Access functions
function get_ae end
"""
    get_aeold(itembuffer::AbstractItemBuffer)

Get the degrees of freedom pertinent to the values object (e.g. cellvalues)
in the itembuffer. 

**Note:** Filled by `NaN`s unless `aold` is passed to `work!`
"""
function get_aeold end

function get_re end
function get_Ke end
function get_material end
function get_values end
"""
    get_old_state(itembuffer::AbstractCellBuffer)

Get the old state variables for the cell. Currently only available for cells and not for facets. 
"""
function get_old_state end

"""
    get_state(itembuffer::AbstractCellBuffer)

Get the state variables for the cell. Currently only available for cells and not for facets. 
"""
function get_state end 

"""
    get_time_increment(itembuffer::AbstractItemBuffer)

Get the time increment to get to the current step. 
Set by [`set_time_increment!`](@ref).
"""
function get_time_increment end 

"""
    get_user_data(itembuffer::AbstractItemBuffer)

Get the `user_data` passed to the `DomainSpec` when setting up the domain. 
This is not modified and always passed around as reference. 
"""
function get_user_data end 

"""
    get_user_cache(itembuffer::AbstractItemBuffer)

Get the `user_cache` created by [`allocate_cell_cache`](@ref).
For multithreaded applications, this cache is copied between 
for each tasks, and can be modified without risking race conditions. 
"""
function get_user_cache end

"""
    get_coupled_buffer(b::AbstractItemBuffer, key::Symbol)

Get the coupled buffer `key` from `b`. To enable this, supply a `coupled_simulations` to
[`work!`](@ref); the coupled itembuffer can be queried just like a normal item buffer,
e.g. by calling `get_state(coupled_buffer)`.
"""
@inline get_coupled_buffer(b::AbstractItemBuffer, key::Symbol) = getfield(get_coupled_buffers(b), key)

get_coupled_buffers(::AbstractItemBuffer) = NamedTuple() # Default: buffer types that don't support coupling

"""
    couple_buffers(itembuffer::AbstractItemBuffer, coupled::Union{CoupledSimulations, NamedTuple})

Refresh `itembuffer`'s coupled-buffer links from `coupled` (either a `CoupledSimulations`, in
sequential work, or a `NamedTuple` of this task's private per-task buffers, in threaded work) and
return `itembuffer`. Called internally, once per `work!` call (before the per-item loop) for
sequential work, or once per task for threaded work. Buffer types that don't support coupling
(e.g. `FacetBuffer`) use this default: a no-op when `coupled` is empty (the common case), or an
error if actually asked to couple.
"""
function couple_buffers(itembuffer::AbstractItemBuffer, coupled)
    _is_empty_coupled(coupled) && return itembuffer
    throw(ArgumentError("$(typeof(itembuffer)) does not support coupled simulations"))
end

_is_empty_coupled(coupled::NamedTuple) = isempty(coupled)
_is_empty_coupled(coupled) = isempty(coupled.sims) # CoupledSimulations

"""
    couple_buffers_or_reuse(itembuffer::AbstractItemBuffer, coupled)

Like [`couple_buffers`](@ref), but skips touching `itembuffer` entirely when there is nothing to
do: `coupled` is empty AND `itembuffer` is already uncoupled. Otherwise (re-)establishes the
links, since `itembuffer` may have been left coupled by a previous, different call.
"""
function couple_buffers_or_reuse(itembuffer::AbstractItemBuffer, coupled)
    (_is_empty_coupled(coupled) && isempty(get_coupled_buffers(itembuffer))) ? itembuffer : couple_buffers(itembuffer, coupled)
end

"""
    Ferrite.celldofs(::AbstractItemBuffer)

Get the degree of freedom indices for the current cell referred to 
by the current item. 
"""
Ferrite.celldofs(::AbstractItemBuffer) = error("Not supported") # Implemented for specific buffers

"""
    Ferrite.getcoordinates(::AbstractItemBuffer)

Get the **cell** coordinates for the current item. 
"""
Ferrite.getcoordinates(::AbstractItemBuffer) = error("Not implemented") # Implemented for specific buffers

"""
    Ferrite.cellid(::AbstractItemBuffer)

Get the cell nr for the current item. 
"""
Ferrite.cellid(::AbstractItemBuffer) = error("Not implemented") # Implemented for specific buffers

"""
    Ferrite.dof_range(::AbstractItemBuffer, ::Symbol)

Get the `dof_range` for a specific field, same as 
`Ferrite.dof_range(::SubDofHandler, ::Symbol)`, 
but fully type-stable. 
"""
Ferrite.dof_range(::AbstractItemBuffer, ::Symbol) = error("Not implemented")

# Set functions
# No docstring - public API is calling this on a domain buffer
function set_time_increment! end

function _replace_material(buf::AbstractItemBuffer, replacement_function)
    new_material = replacement_function(get_material(buf))
    return _replace_material_with(buf, new_material)
end

# Common parts for TaskLocals interface
function scatter!(task::AbstractItemBuffer, base::AbstractItemBuffer)
    set_time_increment!(task, get_time_increment(base))
end
gather!(::AbstractItemBuffer, ::AbstractItemBuffer) = nothing
