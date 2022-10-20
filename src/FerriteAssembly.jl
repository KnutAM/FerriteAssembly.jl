module FerriteAssembly
using Ferrite, ForwardDiff

include("utils.jl")
include("ferrite_additions.jl")
include("CellBuffer.jl")
include("scaling.jl")
include("assembly.jl")
include("states.jl")

export doassemble!, CellBuffer
export create_threaded_CellBuffers, create_threaded_assemblers
export create_states

"""
    element_routine!(
        Ke::AbstractMatrix, re::AbstractVector, 
        ae::AbstractVector, ae_old::AbstractVector,
        state, material, cellvalues, 
        dh_fh::Union{DofHandler,FieldHandler}, Δt, materialcache
        )

The main function to be overloaded for the specific `material` and `cellvalues`.
This function should modify the residual `r` and the element stiffness matrix `Ke`,
such that `Ke=∂(re)/∂(ae)`, where `ae` are the element degrees of freedom. `ae_old` 
are the corresponding old values 
* `state` should contain a state description for the element. 
  Typically, `state` will be a vector with a state variable for each 
  integration point, but it can also be any other type for each element. 
  On input, these are the old values and should be mutated to the updated 
  value for the current time step. 
* The user defined `material` variable usually contain the material parameters
* `cellvalues` should contain the cell values for the given element. 
  It can also be a tuple or named tuple of cellvalues. 
* When the regular `DofHandler` is used, `dh_fh::DofHandler` is passed to the element 
  routine. However, if the `MixedDofHandler` is used, one of its fieldhandlers are passed 
  as `dh_fh::FieldHandler`. This gives the option to call `dof_range(dh_fh, field::Symbol)` 
  for multi-field problems. 
* Finally, a time increment `Δt` and an optional `materialcache` are passed 
  into the element as well.
"""
function element_routine!(
    Ke::AbstractMatrix, re::AbstractVector, 
    ae::AbstractVector, ae_old::AbstractVector, 
    state, material, cellvalues, 
    dh_fh::Union{DofHandler,FieldHandler}, Δt, materialcache
    )
    
    rf!(re_, ae_) = element_residual!(
        re_, ae_, ae_old, state, material, 
        cellvalues, dh_fh, Δt, materialcache
        )
    try
        ForwardDiff.jacobian!(Ke, rf!, re, ae)
    catch e
        if isa(e, MethodError)
            if e.f === element_residual!
                println("Tried to do automatic differentiation to get Ke, but could")
                println("not find a correctly defined `element_residual!` method")
            elseif e.f === convert    
                println("If you get a conversion error, this is likely because")
                println("the element routine is called with dual number inputs.")
                println("If you try to mutate values (e.g. state variables),")
                println("you must use ForwardDiff.value()")
                println("NOTE: If you do this on some values that later affects")
                println("      the output re, the stiffness will be wrong!")
            end
        end
        rethrow(e)
    end
end

"""
    element_residual!(
        re::AbstractVector, 
        ae::AbstractVector, ae_old::AbstractVector,
        state, material, cellvalues, 
        dh_fh::Union{DofHandler,FieldHandler}, Δt, materialcache
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

end
