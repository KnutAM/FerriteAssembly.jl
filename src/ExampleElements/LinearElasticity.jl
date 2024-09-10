@doc raw"""
    ElasticPlaneStrain(;E=2.e3, ν=0.3)

For solving linear elasticity for plane strain, where Young's modulus, `E`, and Poisson's ratio, `ν`,
is used to construct the correct stiffness tensor, ``\boldsymbol{\mathsf{C}}``, 
such that the stress, ``\boldsymbol{\sigma}=\boldsymbol{\mathsf{C}}:\boldsymbol{\epsilon}``, 
where the strain tensor, ``\boldsymbol{\epsilon}=[\boldsymbol{u}\otimes\nabla]^\mathrm{sym}``,
is calculated from the displacement field, ``\boldsymbol{u}(\boldsymbol{x},t)``. 

The strong form of the mechanical quasi-static equilibrium is
```math
    \boldsymbol{\sigma} \cdot \nabla + \boldsymbol{b} = 0 \quad \textbf{x} \in \Omega,
```
with the corresponding weak form, 
```math
   \int_\Omega [\boldsymbol{\delta u}\otimes\nabla]^\mathrm{sym} : \boldsymbol{\sigma} \mathrm{d}\Omega 
   = \int_\Gamma \boldsymbol{\delta u} \cdot \boldsymbol{t} \mathrm{d}\Gamma 
   + \int_\Omega \boldsymbol{\delta u} \cdot \boldsymbol{b} \mathrm{d}\Omega
```
The external loading on the right hand side is not included in the element, but can be implemented 
using `FerriteNeumann.jl`.  (Note that this constructor returns `LinearElastic`, which just stores the 
correct stiffness tensor for the case of isotropic plane strain)
"""
ElasticPlaneStrain(;kwargs...) = LinearElastic(Val(:planestrain); kwargs...)

struct LinearElastic{Dim,T,N}
    C::SymmetricTensor{4,Dim,T,N}
end
function LinearElastic(::Val{Dim}=Val(3);E=2.e3, ν=0.3) where Dim
    G = E / 2(1 + ν)
    K = E / 3(1 - 2ν)
    I2 = one(SymmetricTensor{2,Dim})
    I4vol = I2⊗I2
    I4dev = one(SymmetricTensor{4,Dim}) - I4vol / 3
    return LinearElastic(2G*I4dev + K*I4vol)
end
LinearElastic(::Val{:planestrain}; kwargs...) = LinearElastic(Val(2); kwargs...)
# Plane stress not implemented
#ElasticPlaneStress(;kwargs...) = LinearElastic(Val(:planestress); kwargs...)
#function LinearElastic(::Val{:planestress}; kwargs...)
#    error("Implement me")
#end
# Type-unstable switch could be made for convenience
# LinearElastic(type::Symbol; kwargs...) = LinearElastic(Val(type); kwargs...)

function FerriteAssembly.element_routine!(Ke, re, state, ae, material::LinearElastic, cv::AbstractCellValues, buffer)
    for q_point in 1:getnquadpoints(cv)
        dΩ = getdetJdV(cv, q_point)
        ϵ = function_symmetric_gradient(cv, q_point, ae)
        σ = material.C ⊡ ϵ
        ## Assemble residual contributions
        for i in 1:getnbasefunctions(cv)
            ∇δNu = symmetric(shape_gradient(cv, q_point, i))
            re[i] += (∇δNu ⊡ σ )*dΩ
            ∇δNu_C = ∇δNu ⊡  material.C
            for j in 1:getnbasefunctions(cv)
                ∇Nu = symmetric(shape_gradient(cv, q_point, j))
                Ke[j,i] += (∇δNu_C ⊡ ∇Nu)*dΩ # Since Ke is symmetric, we calculate Ke' to index faster
            end
        end
    end
end

function FerriteAssembly.element_residual!(re, state, ae, material::LinearElastic, cv::AbstractCellValues, buffer)
    for q_point in 1:getnquadpoints(cv)
        dΩ = getdetJdV(cv, q_point)
        ϵ = function_symmetric_gradient(cv, q_point, ae)
        σ = material.C ⊡ ϵ
        ## Assemble residual contributions
        for i in 1:getnbasefunctions(cv)
            ∇δNu = symmetric(shape_gradient(cv, q_point, i))
            re[i] += (∇δNu ⊡ σ )*dΩ
        end
    end
end