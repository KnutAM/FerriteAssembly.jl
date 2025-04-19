#=
# Phase-field fracture

TODO: Describe theory, use micromorphic variable from Ritu...
=#
using Ferrite, FerriteAssembly

@kwdef struct PhaseFieldFracture{C, T}
    G::T    # Elastic shear modulus
    K::T    # Elastic bulk modulus
end
PhaseFieldFracture{C}(;kwargs...) where {C} = PhaseFieldFracture{C, Float64}(;kwargs...)