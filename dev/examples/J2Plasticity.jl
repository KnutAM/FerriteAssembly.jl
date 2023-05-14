# ## Material model
# We use the simple plasticity model from `Ferrite.jl`'s 
# [example](https://ferrite-fem.github.io/Ferrite.jl/stable/examples/plasticity/).
# First, we define the required material struct
struct J2Plasticity{T, S <: SymmetricTensor{4, 3, T}} <: AbstractMaterial
    G::T  # Shear modulus
    K::T  # Bulk modulus
    σ₀::T # Initial yield limit
    H::T  # Hardening modulus
    Dᵉ::S # Elastic stiffness tensor
end
function elastic_conversion(E, ν)
    δ(i,j) = i == j ? 1.0 : 0.0 # helper function
    G = E / 2(1 + ν)
    K = E / 3(1 - 2ν)

    Isymdev(i,j,k,l) = 0.5*(δ(i,k)*δ(j,l) + δ(i,l)*δ(j,k)) - 1.0/3.0*δ(i,j)*δ(k,l)
    temp(i,j,k,l) = 2.0G *( 0.5*(δ(i,k)*δ(j,l) + δ(i,l)*δ(j,k)) + ν/(1.0-2.0ν)*δ(i,j)*δ(k,l))
    Dᵉ = SymmetricTensor{4, 3}(temp)
    return Dᵉ, G, K
end
function J2Plasticity(E, ν, σ₀, H)
    Dᵉ, G, K = elastic_conversion(E, ν)
    return J2Plasticity(G, K, σ₀, H, Dᵉ)
end;
# And the state variable struct 
struct MaterialState{T, S <: SecondOrderTensor{3, T}} <: AbstractMaterialState
    ϵᵖ::S # plastic strain
    σ::S # stress
    k::T # hardening variable
end;
# Then the initial material state function
function MaterialModelsBase.initial_material_state(::J2Plasticity)
    return MaterialState(zero(SymmetricTensor{2,3}), zero(SymmetricTensor{2,3}), 0.0)
end;
# And finally the `material_response` function 
function MaterialModelsBase.material_response(
    material::J2Plasticity, ϵ::SymmetricTensor{2,3}, state::MaterialState, 
    Δt, cache=get_cache(material), args...; kwargs...)
    ## unpack some material parameters
    G = material.G
    H = material.H

    ## We use (•)ᵗ to denote *trial*-values
    σᵗ = material.Dᵉ ⊡ (ϵ - state.ϵᵖ) # trial-stress
    sᵗ = dev(σᵗ)         # deviatoric part of trial-stress
    J₂ = 0.5 * sᵗ ⊡ sᵗ  # second invariant of sᵗ
    σᵗₑ = sqrt(3.0*J₂)   # effective trial-stress (von Mises stress)
    σʸ = material.σ₀ + H * state.k # Previous yield limit

    φᵗ  = σᵗₑ - σʸ # Trial-value of the yield surface

    if φᵗ < 0.0 # elastic loading
        return σᵗ, material.Dᵉ, MaterialState(state.ϵᵖ, σᵗ, state.k)
    else # plastic loading
        h = H + 3G
        μ =  φᵗ / h   # plastic multiplier

        c1 = 1 - 3G * μ / σᵗₑ
        s = c1 * sᵗ           # updated deviatoric stress
        σ = s + vol(σᵗ)       # updated stress

        ## Compute algorithmic tangent stiffness ``D = \frac{\Delta \sigma }{\Delta \epsilon}``
        κ = H * (state.k + μ) # drag stress
        σₑ = material.σ₀ + κ  # updated yield surface

        δ(i,j) = i == j ? 1.0 : 0.0
        Isymdev(i,j,k,l)  = 0.5*(δ(i,k)*δ(j,l) + δ(i,l)*δ(j,k)) - 1.0/3.0*δ(i,j)*δ(k,l)
        Q(i,j,k,l) = Isymdev(i,j,k,l) - 3.0 / (2.0*σₑ^2) * s[i,j]*s[k,l]
        b = (3G*μ/σₑ) / (1.0 + 3G*μ/σₑ)

        Dtemp(i,j,k,l) = -2G*b * Q(i,j,k,l) - 9G^2 / (h*σₑ^2) * s[i,j]*s[k,l]
        D = material.Dᵉ + SymmetricTensor{4, 3}(Dtemp)

        ## Return new state
        Δϵᵖ = 3/2 * μ / σₑ * s # plastic strain
        ϵᵖ = state.ϵᵖ + Δϵᵖ    # plastic strain
        k = state.k + μ        # hardening variable
        return σ, D, MaterialState(ϵᵖ, σ, k)
    end
end;

# Elastic material 
# Then, we also define an elastic material 
struct ElasticMaterial{T<:SymmetricTensor{4,3}} <: AbstractMaterial
    D::T
end
function ElasticMaterial(;E, ν)
    D, _, _ = elastic_conversion(E, ν)
    return ElasticMaterial(D)
end;

# This material requires not state, so `MaterialModelsBase.jl`'s `NoMaterialState`
# will be created by default. Hence, we only need to define the `material_response`
function MaterialModelsBase.material_response(
    material::ElasticMaterial, ϵ::SymmetricTensor{2,3}, state, 
    Δt, cache=get_cache(material), args...; kwargs...)
    σ = material.D ⊡ ϵ
    return σ, material.D, NoMaterialState()
end;