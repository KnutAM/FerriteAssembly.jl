@doc raw"""
    PoroElasticPlaneStrain(;E=2.e3, ν=0.3, k=0.05, α=1.0, β=1/2e3)

The strong forms are given as
```math
\begin{aligned}
\boldsymbol{\sigma}(\boldsymbol{\epsilon}, p) \cdot \boldsymbol{\nabla} &= \boldsymbol{0} \\
\dot{\Phi}(\boldsymbol{\epsilon}, p) + \boldsymbol{w}(p) \cdot \boldsymbol{\nabla} &= 0
\end{aligned}
```
where 
``\boldsymbol{\epsilon} = \left[\boldsymbol{u}\otimes\boldsymbol{\nabla}\right]^\mathrm{sym}`` 
The constitutive relationships are 
```math
\begin{aligned}
\boldsymbol{\sigma} &= \boldsymbol{\mathsf{E}}:\boldsymbol{\epsilon} - \alpha p \boldsymbol{I} \\
\boldsymbol{w} &= - k \boldsymbol{\nabla} p \\
\Phi &= \phi + \alpha \mathrm{tr}(\boldsymbol{\epsilon}) + \beta p
\end{aligned}
``` 
with 
``\boldsymbol{\mathsf{E}}=2G \boldsymbol{\mathsf{I}}^\mathrm{dev} + 3K \boldsymbol{I}\otimes\boldsymbol{I}``.
The material parameters are then the 
shear modulus, ``G``, 
bulk modulus, ``K``, 
permeability, ``k``,  
Biot's coefficient, ``\alpha``, and
liquid compressibility, ``\beta``.
The porosity, ``\phi``, doesn't enter into the equations 
(A different porosity leads to different skeleton stiffness and permeability).

The weak forms are
```math
\begin{aligned}
\int_\Omega \left[\left[\boldsymbol{\delta u}\otimes\boldsymbol{\nabla}\right]^\mathrm{sym}:
\boldsymbol{\mathsf{E}}:\boldsymbol{\epsilon} - \boldsymbol{\delta u} \cdot \boldsymbol{\nabla} \alpha p\right] \mathrm{d}\Omega 
&= \int_\Gamma \boldsymbol{\delta u} \cdot \boldsymbol{t} \mathrm{d} \Gamma \\
\int_\Omega \left[\delta p \left[\alpha \dot{\boldsymbol{u}} \cdot \boldsymbol{\nabla} + \beta \dot{p}\right] + 
\boldsymbol{\nabla}(\delta p) \cdot [k \boldsymbol{\nabla}]\right] \mathrm{d}\Omega 
&= -\int_\Gamma \delta p w_\mathrm{n} \mathrm{d} \Gamma 
\end{aligned}
```
where ``\boldsymbol{t}=\boldsymbol{n}\cdot\boldsymbol{\sigma}`` is the traction and 
``w_\mathrm{n} = \boldsymbol{n}\cdot\boldsymbol{w}`` is the normal flux.  
"""
struct PoroElasticPlaneStrain{T}
    C::SymmetricTensor{4,2,T,9} # 2D plane strain elastic stiffness
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
    C = 2G*I4dev + K*I4vol
    return PoroElasticPlaneStrain(C, k, α, β)
end

function FerriteAssembly.element_residual!(re, new_state, ae, material::PoroElasticPlaneStrain, cv::NamedTuple, buffer)
    ae_old = FerriteAssembly.get_aeold(buffer)
    Δt = FerriteAssembly.get_time_increment(buffer)
    udofs = dof_range(buffer, :u)
    pdofs = dof_range(buffer, :p)

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
        σeff = material.C ⊡ ϵ

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