abstract type AbstractItemBuffer end

# Access functions
"""
    get_dofbuffer(itembuffer::AbstractItemBuffer)

Get the [`DofBuffer`](@ref FerriteAssembly.DofBuffer) for the item. Internal function, which each
item buffer must define. The accessors for the local dofs, their values, and the local residual
and stiffness are defined for all item buffers by using this function.
"""
function get_dofbuffer end

"""
    get_ae(itembuffer::AbstractItemBuffer)

Get the values of the item's local degrees of freedom stored in the itembuffer,
e.g. for a coupled buffer, see [`get_coupled_buffer`](@ref).

**Note:** In `element_residual!`, use the input `ae` instead, as `get_ae(itembuffer)`
doesn't contain the dual numbers required for automatic differentiation.
"""
@inline get_ae(b::AbstractItemBuffer) = get_dofbuffer(b).ae
"""
    get_aeold(itembuffer::AbstractItemBuffer)

Get the degrees of freedom pertinent to the values object (e.g. cellvalues)
in the itembuffer. 

**Note:** Filled by `NaN`s unless `aold` is passed to `work!`
"""
@inline get_aeold(b::AbstractItemBuffer) = get_dofbuffer(b).ae_old

@inline get_re(b::AbstractItemBuffer) = get_dofbuffer(b).re
@inline get_Ke(b::AbstractItemBuffer) = get_dofbuffer(b).Ke
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

Get the coupled buffer `key` from `b`. To enable this, build a [`CoupledSimulations`](@ref)
group and `work!` its member handles. The coupled buffer can be queried just like a normal
item buffer, e.g. by calling `get_state(coupled_buffer)`.
"""
@inline get_coupled_buffer(b::AbstractItemBuffer, key::Symbol) = getfield(get_coupled_buffers(b), key)

"""
    Ferrite.celldofs(::AbstractItemBuffer)

Get the degree of freedom indices for the current cell referred to 
by the current item. For a domain with algebraic variables 
(see the `algebraic_per_cell` keyword of [`DomainSpec`](@ref)), these 
are followed by the degree of freedom indices for the algebraic variables.
"""
@inline Ferrite.celldofs(b::AbstractItemBuffer) = get_dofbuffer(b).dofs

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
but fully type-stable. For a domain with algebraic variables 
(see the `algebraic_per_cell` keyword of [`DomainSpec`](@ref)),
the name of an algebraic variable gives the local dof range for that variable.
"""
@inline Ferrite.dof_range(b::AbstractItemBuffer, name::Symbol) = get_dofbuffer(b).dofrange[name]

# Internal overload for now
Ferrite.getfieldnames(b::AbstractItemBuffer) = keys(get_dofbuffer(b).dofrange)

# The `AlgebraicLayout` for the buffer, or `nothing` if no algebraic variables are used
get_layout(::AbstractItemBuffer) = nothing

has_algebraic_dofs(b::AbstractItemBuffer) = get_layout(b) !== nothing

"""
    get_algebraic_values(itembuffer::AbstractItemBuffer, name::Symbol)

Get the `Ferrite.AlgebraicValues` for the algebraic variable `name`, which must be one of the
variables given as `algebraic_per_cell` to the [`DomainSpec`](@ref). Together with
`dof_range(itembuffer, name)`, the value of the variable is then
`algebraic_value(get_algebraic_values(itembuffer, name), ae, dof_range(itembuffer, name))`.
"""
@inline get_algebraic_values(b::AbstractItemBuffer, name::Symbol) = _get_algebraic_values(get_layout(b), name)

# The names of the spatial fields (i.e. excluding the algebraic variables in `getfieldnames`)
_spatial_fieldnames(b::AbstractItemBuffer) = _spatial_fieldnames(Ferrite.getfieldnames(b), get_layout(b))

# Write the current global dofs of the algebraic variables to the local dofs of the buffer
update_algebraic_dofs!(b::AbstractItemBuffer, dh) = update_algebraic_dofs!(b, get_layout(b), dh)
update_algebraic_dofs!(::AbstractItemBuffer, ::Nothing, dh) = nothing
update_algebraic_dofs!(b::AbstractItemBuffer, layout::AlgebraicLayout, dh) = update_algebraic_dofs!(celldofs(b), layout, dh)

# Copy the global dofs of the algebraic variables from the local dofs of `src` to those of `dst`
copy_algebraic_dofs!(dst::AbstractItemBuffer, src::AbstractItemBuffer) = copy_algebraic_dofs!(dst, src, get_layout(src))
copy_algebraic_dofs!(::AbstractItemBuffer, ::AbstractItemBuffer, ::Nothing) = nothing
function copy_algebraic_dofs!(dst::AbstractItemBuffer, src::AbstractItemBuffer, layout::AlgebraicLayout)
    return copy_algebraic_dofs!(celldofs(dst), celldofs(src), layout)
end

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
    copy_algebraic_dofs!(task, base)
end
gather!(::AbstractItemBuffer, ::AbstractItemBuffer) = nothing
