# # `MaterialModelsBase.jl`: Plasticity
# This example shows how any material following the 
# [`MaterialModelsBase.jl`](https://github.com/KnutAM/MaterialModelsBase.jl)
# interface can be assembled with `FerriteAssembly.jl`. 
# We start by the required packages
using Tensors, MaterialModelsBase, Ferrite, FerriteAssembly

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
function J2Plasticity(E, ν, σ₀, H)
    δ(i,j) = i == j ? 1.0 : 0.0 # helper function
    G = E / 2(1 + ν)
    K = E / 3(1 - 2ν)

    Isymdev(i,j,k,l) = 0.5*(δ(i,k)*δ(j,l) + δ(i,l)*δ(j,k)) - 1.0/3.0*δ(i,j)*δ(k,l)
    temp(i,j,k,l) = 2.0G *( 0.5*(δ(i,k)*δ(j,l) + δ(i,l)*δ(j,k)) + ν/(1.0-2.0ν)*δ(i,j)*δ(k,l))
    Dᵉ = SymmetricTensor{4, 3}(temp)
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

# ## Element routine
# The element routine for any material that follows the `MaterialModelsBase.jl` interface 
# can be coded as. The boundary traction is not included, as this can be handled separately
# with [FerriteNeumann.jl](https://github.com/KnutAM/FerriteNeumann.jl)
function FerriteAssembly.element_routine!(
    Ke::AbstractMatrix, re::AbstractVector, state::Vector{<:AbstractMaterialState}, 
    ae::AbstractVector, material::AbstractMaterial, cellvalues::CellVectorValues, dh_fh, Δt, buffer::CellBuffer)
    cache = buffer.cache
    n_basefuncs = getnbasefunctions(cellvalues)
    for q_point in 1:getnquadpoints(cellvalues)
        ## For each integration point, compute stress and material stiffness
        ϵ = function_symmetric_gradient(cellvalues, q_point, ae) # Total strain
        σ, D, state[q_point] = material_response(material, ϵ, state[q_point], Δt, cache)

        dΩ = getdetJdV(cellvalues, q_point)
        for i in 1:n_basefuncs
            δϵ = shape_symmetric_gradient(cellvalues, q_point, i)
            re[i] += (δϵ ⊡ σ) * dΩ ## add internal force to residual
            for j in 1:n_basefuncs
                Δϵ = shape_symmetric_gradient(cellvalues, q_point, j)
                Ke[i, j] += δϵ ⊡ D ⊡ Δϵ * dΩ
            end
        end
    end
end;

# ## Assembly
# With all required functions defined, we can now setup and assemble the finite element problem 
material = J2Plasticity(200.0e9, 0.3, 200.0e6, 10.0e9);
grid = generate_grid(Tetrahedron, (20,2,4), zero(Vec{3}), Vec((10.0,1.0,1.0)));
cellvalues = CellVectorValues(
    QuadratureRule{3,RefTetrahedron}(2), Lagrange{3, RefTetrahedron, 1}());
dh = DofHandler(grid); push!(dh, :u, 3); close!(dh); # Create dofhandler
K = create_sparsity_pattern(dh);
r = zeros(ndofs(dh));

# Using the `create_states` function, we can easily create the storage of state variables 
# suitable for the chosen dof handler 
states = create_states(dh, _->initial_material_state(material), cellvalues);

# And then we create the buffer for saving cell-related variables
buffer = CellBuffer(dh::DofHandler, cellvalues, material, nothing, get_cache(material));

# And for an initial guess of displacements, `a`, 
# then create the assembler and do the assembly
a = zeros(ndofs(dh))
assembler = start_assemble(K,r)
doassemble!(assembler, buffer, states, dh, a);