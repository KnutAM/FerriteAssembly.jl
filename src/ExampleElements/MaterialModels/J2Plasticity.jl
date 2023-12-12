# Material struct and parameters
"""
    J2Plasticity(;E, ν, σ0, H) <: MaterialModelsBase.AbstractMaterial

This plasticity model is taken from `Ferrite.jl`'s plasticity 
[example](https://ferrite-fem.github.io/Ferrite.jl/v0.3.14/examples/plasticity/),
and considers linear isotropic hardening, with Young's modulus, `E`, Poisson's
ratio, `ν`, initial yield limit, `σ0`, and hardening modulus, `H`. It is defined 
as an `AbstractMaterial` following the `MaterialModelsBase` interface. 
"""
struct J2Plasticity{T, S <: SymmetricTensor{4, 3, T}} <: AbstractMaterial
    G::T  # Shear modulus
    K::T  # Bulk modulus
    σ0::T # Initial yield limit
    H::T  # Hardening modulus
    D::S  # Elastic stiffness tensor
end
function J2Plasticity(;E, ν, σ0, H)
    δ(i,j) = i == j ? 1.0 : 0.0 # helper function
    G = E / 2(1 + ν)
    K = E / 3(1 - 2ν)

    Isymdev(i,j,k,l) = 0.5*(δ(i,k)*δ(j,l) + δ(i,l)*δ(j,k)) - 1.0/3.0*δ(i,j)*δ(k,l)
    temp(i,j,k,l) = 2.0G *( 0.5*(δ(i,k)*δ(j,l) + δ(i,l)*δ(j,k)) + ν/(1.0-2.0ν)*δ(i,j)*δ(k,l))
    D = SymmetricTensor{4, 3}(temp)
    return J2Plasticity(G, K, σ0, H, D)
end

# State variables
struct J2PlasticityState{T, S <: SecondOrderTensor{3, T}} <: AbstractMaterialState
    ϵp::S # plastic strain
    k::T # hardening variable
end
function MaterialModelsBase.initial_material_state(::J2Plasticity)
    return J2PlasticityState(zero(SymmetricTensor{2,3}), 0.0)
end

# The main `material_response` function 
function MaterialModelsBase.material_response(
    material::J2Plasticity, ϵ::SymmetricTensor{2,3}, state::J2PlasticityState, Δt, cache, args...)
    ## unpack some material parameters
    G = material.G
    H = material.H

    ## We use (•)ᵗ to denote *trial*-values
    σᵗ = material.D ⊡ (ϵ - state.ϵp) # trial-stress
    sᵗ = dev(σᵗ)         # deviatoric part of trial-stress
    J₂ = 0.5 * sᵗ ⊡ sᵗ  # second invariant of sᵗ
    σᵗₑ = sqrt(3.0*J₂)   # effective trial-stress (von Mises stress)
    σʸ = material.σ0 + H * state.k # Previous yield limit

    φᵗ  = σᵗₑ - σʸ # Trial-value of the yield surface

    if φᵗ < 0.0 # elastic loading
        return σᵗ, material.D, state
    else # plastic loading
        h = H + 3G
        μ =  φᵗ / h   # plastic multiplier

        c1 = 1 - 3G * μ / σᵗₑ
        s = c1 * sᵗ           # updated deviatoric stress
        σ = s + vol(σᵗ)       # updated stress

        ## Compute algorithmic tangent stiffness ``D = \frac{\Delta \sigma }{\Delta \epsilon}``
        κ = H * (state.k + μ) # drag stress
        σₑ = material.σ0 + κ  # updated yield surface

        δ(i,j) = i == j ? 1.0 : 0.0
        Isymdev(i,j,k,l)  = 0.5*(δ(i,k)*δ(j,l) + δ(i,l)*δ(j,k)) - 1.0/3.0*δ(i,j)*δ(k,l)
        Q(i,j,k,l) = Isymdev(i,j,k,l) - 3.0 / (2.0*σₑ^2) * s[i,j]*s[k,l]
        b = (3G*μ/σₑ) / (1.0 + 3G*μ/σₑ)

        Dtemp(i,j,k,l) = -2G*b * Q(i,j,k,l) - 9G^2 / (h*σₑ^2) * s[i,j]*s[k,l]
        D = material.D + SymmetricTensor{4, 3}(Dtemp)

        ## Return new state
        Δϵᵖ = 3/2 * μ / σₑ * s # plastic strain
        ϵp = state.ϵp + Δϵᵖ    # plastic strain
        k = state.k + μ        # hardening variable
        return σ, D, J2PlasticityState(ϵp, k)
    end
end
