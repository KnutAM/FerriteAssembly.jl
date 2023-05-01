# ### Poroelastic material
# For the poroelastic material, we reuse the elastic part from above, 
# but add additionally required properties 
struct PoroElasticPlaneStrain{T}
    dσdϵ::SymmetricTensor{4,2,T,9} # 2D plane strain elastic stiffness
    k::T    # [mm^4/Ns] Permeability
    α::T    # [-] Biot's coefficient
    β::T    # [1/MPa] Liquid bulk modulus
end 
function PoroElasticPlaneStrain(;E=2.e3, ν=0.3, k=0.05, α=1.0, β=1/2e3)
    G = E / 2(1 + ν)
    K = E / 3(1 - 2ν)
    I2 = one(SymmetricTensor{2,2})
    I4vol = I2⊗I2
    I4dev = minorsymmetric(otimesu(I2,I2)) - I4vol / 3
    dσdϵ = 2G*I4dev + K*I4vol
    return PoroElasticPlaneStrain(dσdϵ, k, α, β)
end

function FerriteAssembly.element_residual!(re, state, ae, material::PoroElasticPlaneStrain, cv::NamedTuple, dh_fh, Δt, buffer)
    ae_old = FerriteAssembly.get_aeold(buffer)
    udofs = dof_range(dh_fh, :u)
    pdofs = dof_range(dh_fh, :p)

    ## Assemble stiffness and force vectors
    for q_point in 1:getnquadpoints(cv[:u])   
        ## Calculate variables in the current quadrature point 
        dΩ = getdetJdV(cv[:u], q_point)
        ϵ = function_symmetric_gradient(cv[:u], q_point, ae, udofs)
        ϵ_old = function_symmetric_gradient(cv[:u], q_point, ae_old, udofs)
        p = function_value(cv[:p], q_point, ae, pdofs)
        p_old = function_value(cv[:p], q_point, ae_old, pdofs)
        ∇p = function_gradient(cv[:p], q_point, ae, pdofs)
        pdot = (p-p_old)/Δt 
        div_udot = (tr(ϵ)-tr(ϵ_old))/Δt
        σeff = material.dσdϵ ⊡ ϵ

        ## Assemble residual contributions
        for (iᵤ, i) in enumerate(udofs)
            ∇δNu = shape_symmetric_gradient(cv[:u], q_point, iᵤ)
            div_δNu = shape_divergence(cv[:u], q_point, iᵤ)
            re[i] += (∇δNu ⊡ σeff - div_δNu*material.α*p)*dΩ
        end
        for (iₚ, i) in enumerate(pdofs)
            δNp = shape_value(cv[:p], q_point, iₚ)
            ∇δNp = shape_gradient(cv[:p], q_point, iₚ)
            re[i] += (δNp*(material.α*div_udot + material.β*pdot) + (∇δNp ⋅ ∇p)*material.k) * dΩ
        end
    end
end;