module FerriteAssembly
using Ferrite, ForwardDiff

include("Multithreading/TaskLocals.jl") # Task-local storage model 
include("Multithreading/TaskChunks.jl") # Thread-safe iteration over chunks of cells

include("Utils/FerriteAddons.jl")       # Temporary solutions until Ferrite is updated
include("Utils/utils.jl")
include("Utils/scaling.jl")
include("states.jl")

include("ItemBuffers/AbstractItemBuffer.jl")

include("DomainBuffers.jl")
include("setup.jl")

include("ItemBuffers/CellBuffer.jl")
include("ItemBuffers/FacetBuffer.jl")
include("Autodiff/autodiff.jl")

include("work.jl")
include("Workers/Assemblers.jl")
include("Workers/Integrators.jl")
include("Workers/QuadratureEvaluator.jl")

include("LoadHandler/LoadHandler.jl")


# Setup 
export DomainSpec, setup_domainbuffer, setup_domainbuffers
# Main functions to use during simulations
export work!, update_states!, set_time_increment!
# Workers
export ReAssembler, KeReAssembler   # Assemblers
export Integrator, SimpleIntegrator # Integrators
# Builtin convenience 
export LoadHandler, Neumann, BodyLoad, DofLoad
export ElementResidualScaling, reset_scaling!

"""
    element_routine!(
        Ke::AbstractMatrix, re::AbstractVector, state,
        ae::AbstractVector, material, cellvalues, buffer)

The main function to be overloaded for the specific `material`. In most cases,
the same implementation can be used for different cellvalues 
(e.g. for different interpolation orders)
This function should modify the element stiffness matrix `Ke`, the residual `re`,
and potentially `state`. The element degree of freedom values, `ae`, are filled by 
`NaN`s unless `a` is passed to [`work!`](@ref).

The following variables can be obtained from `buffer`.
* [`get_old_state(buffer)`](@ref FerriteAssembly.get_old_state)
* [`get_aeold(buffer)`](@ref FerriteAssembly.get_aeold)
* [`get_time_increment(buffer)`](@ref FerriteAssembly.get_time_increment)
* [`dof_range(buffer, fieldname::Symbol)`](@ref Ferrite.dof_range)
* [`getcoordinates(buffer)`](@ref Ferrite.getcoordinates)
* [`celldofs(buffer)`](@ref Ferrite.celldofs)
* [`cellid(buffer)`](@ref Ferrite.cellid)
* [`get_user_data(buffer)`](@ref FerriteAssembly.get_user_data)
* [`get_user_cache(buffer)`](@ref FerriteAssembly.get_user_cache)
"""
element_routine!(args...) = element_routine_ad!(args...)    # If not defined, try to use automatic differentiation

"""
    element_residual!(
        re::AbstractVector, state,
        ae::AbstractVector, material, cellvalues, buffer)

To calculate the element tangent stiffness `Ke` automatically by using `ForwardDiff`,
it is possible to overload `element_residual!` instead of `element_routine!`. See 
[`element_routine!`](@ref) for a description of the input parameters. 

!!! note 
    `MethodError` with `ForwardDiff.Dual`\\
    When using automatic differentiation for elements with state variables (or other mutating values in e.g. cache),
    an error will be thrown if trying to change the type in many cases. 
    When mutating `state`, call `ForwardDiff.value()` on the value to be assigned **after** 
    it will no longer be used to calculate `re`. (If done on a value which is later affects `re` inside the element,
    the tangent, `Ke`, will be wrong.)
"""
function element_residual! end

"""
    facet_routine!(Ke, re, ae, material, facetvalues, facetbuffer)

Calculate contributions to the stiffness matrix and residual vector from a 
facet domain. It can be used, for example, to implement Robin boundary conditions.
"""
function facet_routine! end

"""
    facet_residual!(re, ae, material, facetvalues, facetbuffer)

Calculate contributions to the residual vector from a facet domain. 
Internally, this is used for calculating Neumann boundary conditions. 
"""
function facet_residual! end 


include("Utils/MaterialModelsBase.jl")
include("ExampleElements/ExampleElements.jl")

end
