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

