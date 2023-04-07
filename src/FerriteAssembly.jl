module FerriteAssembly
using Ferrite, ForwardDiff

include("utils.jl")
include("ferrite_additions.jl")
include("scaling.jl")
include("CellBuffer.jl")
include("autodiff.jl")
include("assembly.jl")
include("states.jl")


export doassemble!
export CellBuffer, AutoDiffCellBuffer, getCellBuffer
export setup_cellbuffer, setup_ad_cellbuffer
export create_states
export ElementResidualScaling, reset_scaling!
export create_threaded_CellBuffers, create_threaded_assemblers, create_threaded_scalings

"""
    element_routine!(
        Ke::AbstractMatrix, re::AbstractVector, state,
        ae::AbstractVector, material, cellvalues, 
        dh_fh, Δt, buffer
        )

The main function to be overloaded for the specific `material` and `cellvalues`.
This function should modify the element stiffness matrix `Ke` and the residual `re`.
* `state` should contain a state description for the element. 
  Typically, `state` will be a vector with a state variable for each 
  integration point, but it can also be any other type for each element. 
  On input, these are the old values and should be mutated to the updated 
  value for the current time step and guess for `ae`. 
* The user defined `material` variable usually contain the material parameters. 
* `cellvalues` should contain the `CellValues` for the given element. 
  It can also be a tuple or named tuple of cellvalues. 
* When the regular `DofHandler` is used, `dh_fh::DofHandler` is passed to the element 
  routine. However, if the `MixedDofHandler` is used, one of its fieldhandlers are passed 
  as `dh_fh::FieldHandler`. This gives the option to call `dof_range(dh_fh, field::Symbol)` 
  for multi-field problems. **Note:** Please do not rely on `dh_fh` for anything but `dof_range`,
  as `dh_fh` may be replaced with another type that only supports `dof_range`.
* `Δt` is time increment given to `doassemble`
* `buffer::Union{CellBuffer, AutoDiffCellBuffer}` can be used to get 
  - `getcoordinates(buffer)::Vector{Vec}`: The cell's coordinates
  - `celldofs(buffer)::Vector{Int}`: The cell's global degrees of freedom numbers
  - `FerriteAssembly.get_aeold(buffer)`: `aold[celldofs]` (if `aold::Nothing` is passed to 
    `doassemble`, a vector, `[NaN for d in celldofs]` is returned)
  - `FerriteAssembly.get_load(buffer)`: The `cell_load` in the buffer, typically used for 
     body loads or source terms. 
  - `FerriteAssembly.get_cache(buffer)`: The `cache` in the buffer - typically used to gather all 
     preallocations if such are necessary
"""
element_routine!(args...) = element_routine_ad!(args...)    # If not defined, try to use automatic differentiation

"""
    element_residual!(
        re::AbstractVector, state, 
        ae::AbstractVector, material, cellvalues, 
        dh_fh, Δt, buffer
        )

To calculate the element tangent stiffness `Ke` automatically by using `ForwardDiff`,
it is possible to overload `element_residual!` instead of `element_routine!`. See 
[`element_routine!`](@ref) for a description of the input parameters. 

Note that in order for this function to work, care must be taken when mutating values
to not change their types. When mutating the `state`, ensure to call `ForwardDiff.value()`
on those values. *Warning*: Only do this at a point when the calculation of `re` is 
unaffected by the values in `state`, otherwise `Ke` will be wrong. 
"""
function element_residual! end


include("ExampleElements/ExampleElements.jl")

end
