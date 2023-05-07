import MaterialModelsBase as MMB

@doc raw"""
    FerriteAssembly.element_routine!(Ke, re, state, ae, m::MaterialModelsBase.AbstractMaterial, args...)

Solve the weak form 
```math
   \int_\Omega [\boldsymbol{\delta u}\otimes\nabla]^\mathrm{sym} : \boldsymbol{\sigma} \mathrm{d}\Omega 
   = \int_\Gamma \boldsymbol{\delta u} \cdot \boldsymbol{t} \mathrm{d}\Gamma 
   + \int_\Omega \boldsymbol{\delta u} \cdot \boldsymbol{b} \mathrm{d}\Omega
```
where ``\sigma`` is calculated with the `material_response` function from 
[`MaterialModelsBase.jl`](https://github.com/KnutAM/MaterialModelsBase.jl). 
Note that `create_cell_state` is already implemented for `<:AbstractMaterial`. 
"""
function FerriteAssembly.element_routine!(
    Ke, re, state_new::Vector{<:MMB.AbstractMaterialState},
    ae, material::MMB.AbstractMaterial, cellvalues::CellVectorValues, buffer)
    cache = FerriteAssembly.get_cache(buffer)
    Δt = FerriteAssembly.get_time_increment(buffer)
    state_old = FerriteAssembly.get_state_old(buffer)
    n_basefuncs = getnbasefunctions(cellvalues)
    for q_point in 1:getnquadpoints(cellvalues)
        # For each integration point, compute stress and material stiffness
        ϵ = function_symmetric_gradient(cellvalues, q_point, ae) # Total strain
        σ, D, state_new[q_point] = MMB.material_response(material, ϵ, state_old[q_point], Δt, cache)

        dΩ = getdetJdV(cellvalues, q_point)
        for i in 1:n_basefuncs
            ∇δN = shape_symmetric_gradient(cellvalues, q_point, i)
            re[i] += (∇δN ⊡ σ) * dΩ # add internal force to residual
            ∇δN_D = ∇δN ⊡ D         # temporary value for speed
            for j in 1:n_basefuncs
                ∇N = shape_symmetric_gradient(cellvalues, q_point, j)
                Ke[i, j] += (∇δN_D ⊡ ∇N) * dΩ
            end
        end
    end
end

function FerriteAssembly.create_cell_state(m::MMB.AbstractMaterial, cv::CellVectorValues, args...)
    return [MMB.initial_material_state(m) for _ in 1:getnquadpoints(cv)]
end