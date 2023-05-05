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
    Ke, re, state::Vector{<:MMB.AbstractMaterialState},
    ae, material::MMB.AbstractMaterial, cellvalues::CellVectorValues, 
    dh_fh, Δt, cb)
    buffer = FerriteAssembly.getCellBuffer(cb)
    cache = buffer.cache
    n_basefuncs = getnbasefunctions(cellvalues)
    for q_point in 1:getnquadpoints(cellvalues)
        # For each integration point, compute stress and material stiffness
        ϵ = function_symmetric_gradient(cellvalues, q_point, ae) # Total strain
        σ, D, state[q_point] = MMB.material_response(material, ϵ, state[q_point], Δt, cache)

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

function FerriteAssembly.element_routine!(Ke, re, state::Nothing, ae, ::MMB.AbstractMaterial, args...; kwargs...)
    msg = "When using MaterialModelsBase materials, state::Vector{<:AbstractMaterialState} is required.\n"*
          "Perhaps you forgot to add material and cellvalues when calling create_states?"
    throw(ArgumentError(msg))
end
