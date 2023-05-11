module FerriteAssembly
using Ferrite, ForwardDiff

# Temporary solutions until Ferrite is updated
include("SubDofHandler.jl")
include("TaskLocals.jl")    # Task-local storage model 
include("utils.jl")

include("scaling.jl")
include("states.jl")
include("CellBuffer.jl")
include("autodiff.jl")

include("Assemblers.jl")
include("setup.jl")
include("assembly.jl")

export setup_assembly, AssemblyDomain
export doassemble!, update_states!

export ElementResidualScaling, reset_scaling!

"""
    element_routine!(
        Ke::AbstractMatrix, re::AbstractVector, state_new,
        ae::AbstractVector, material, cellvalues, buffer)

The main function to be overloaded for the specific `material`. In most cases,
the same implementation can be used for different cellvalues 
(e.g. for different interpolation orders)
This function should modify the element stiffness matrix `Ke`, the residual `re`,
and potentially `state_new`. The element degree of freedom values, `ae`, are filled by 
`NaN`s unless `a` is passed to [`doassemble!`](@ref).

The following variables can be obtained from `buffer`.
* [`get_state_old(buffer)`](@ref FerriteAssembly.get_state_old)
* [`get_aeold(buffer)`](@ref FerriteAssembly.get_aeold)
* [`get_time_increment(buffer)`](@ref FerriteAssembly.get_time_increment)
* [`dof_range(buffer, fieldname::Symbol)`](@ref Ferrite.dof_range)
* [`getcoordinates(buffer)`](@ref Ferrite.getcoordinates)
* [`celldofs(buffer)`](@ref Ferrite.celldofs)
* [`cellid(buffer)`](@ref Ferrite.cellid)
* [`get_user_data(buffer)`](@ref FerriteAssembly.get_user_data)
* [`get_cache(buffer)`](@ref FerriteAssembly.get_cache)
"""
element_routine!(args...) = element_routine_ad!(args...)    # If not defined, try to use automatic differentiation

"""
    element_residual!(
        re::AbstractVector, state_new,
        ae::AbstractVector, material, cellvalues, buffer)

To calculate the element tangent stiffness `Ke` automatically by using `ForwardDiff`,
it is possible to overload `element_residual!` instead of `element_routine!`. See 
[`element_routine!`](@ref) for a description of the input parameters. 

!!! note `MethodError` with `ForwardDiff.Dual`
    When using automatic differentiation for elements with state variables (or other mutating values in e.g. cache),
    an error will be thrown if trying to change the type in many cases. 
    When mutating `new_state`, call `ForwardDiff.value()` on the value to be assigned **after** 
    it will no longer be used to calculate `re`. (If done on a value which is later affects `re` inside the element,
    the tangent, `Ke`, will be wrong.)
"""
function element_residual! end

include("MaterialModelsBase.jl")
include("ExampleElements/ExampleElements.jl")

end
