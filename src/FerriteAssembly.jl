module FerriteAssembly
using Ferrite, ForwardDiff

include("CellCache.jl")
include("scaling.jl")
include("assembly.jl")

export doassemble!, CellCache

function element_routine!(
    Ke::AbstractMatrix, re::AbstractVector, 
    ae_new::AbstractVector, ae_old::AbstractVector,
    state::AbstractVector, material, cellvalues, 
    dh_fh::Union{DofHandler,FieldHandler}, Δt, materialcache
    )
    
    rf!(re_, ae_) = element_residual!(
        re_, ae_, ae_old, state, material, 
        cellvalues, dh_fh, Δt, materialcache
        )
    try
        ForwardDiff.jacobian!(Ke, rf!, re, ae_new)
    catch e
        if isa(e, MethodError)    
            println("If you get a conversion error, this is likely because")
            println("the element routine is called with dual number inputs.")
            println("If you try to mutate values (e.g. state variables),")
            println("you must use ForwardDiff.value()")
            println("NOTE: If you do this on some values that later affects")
            println("      the output re, the stiffness will be wrong!")
        end
        rethrow(e)
    end
end

function element_residual!(
    re::AbstractVector, ae_new::AbstractVector, ae_old::AbstractVector,
    state::AbstractVector, material, cellvalues, 
    dh_fh::Union{DofHandler,FieldHandler}, Δt, materialcache
    )
    println("Only the interface function definition is being called!")
    throw(MethodError(
        element_residual!, 
        (re, ae_new, ae_old, state, material, cellvalues, dh_fh, Δt, materialcache)
        ))
end

end
