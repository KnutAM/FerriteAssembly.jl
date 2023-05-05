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
        σ, D, state[q_point] = material_response(material, ϵ, state[q_point], Δt, cache)

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


# =================================================== #
@doc raw"""
    WeakForm(f::Function)

Solve the problem with primary variable, ``u``, variation, ``\delta u``, and a weak form
```math
   \int_\Omega f(\delta u, \delta u \otimes\nabla, u, u \otimes\nabla, \dot{u}, \dot{u} \otimes\nabla) \mathrm{d}\Omega 
   - \int_\Gamma \delta u\ h(x,t,n)\ \mathrm{d}\Gamma - \int_\Omega \delta u\ b(x,t)\ \mathrm{d}\Omega = 0
```
where the function `f` is given to the weak form, and `h` and `b` are given with FerriteNeumann.
This is only intended for testing, and is not optimized for speed 
(`f` is called for each shape function for each integration point)

# Examples 
## Transient heat flow
*Weak form*
```math
    \int_\Omega \delta u\ c\ \dot{u} - k\ [\nabla \delta u] \cdot [\nabla u]\ \mathrm{d}\Omega 
   - \int_\Gamma \delta u\ q_\mathrm{n}\ \mathrm{d}\Gamma + \int_\Omega \delta u\ b\ \mathrm{d}\Omega = 0
```
*Implementation*
```julia
c = 1.0; k = 1.0; # heat capacity and heat conductivity (material parameters)
qn = 1.0; b=1.0;  # Normal boundary flux and internal heat source (external loading)
material = WeakForm((δu, ∇δu, u, ∇u, u_dot, ∇u_dot) -> δu*c*u_dot + k*(∇δu ⋅ ∇u))
nh = NeumannHandler(dh)
add!(nh, Neumann(:u, 2, getfaceset(dh.grid, "right"), (x,t,n)->qn))
add!(nh, BodyLoad(:c, 1, (x,t)->b))
```

## Linear elasticity
*Weak form*
```math
   \int_\Omega [\boldsymbol{\delta u}\otimes\nabla]^\mathrm{sym} : \boldsymbol{\sigma}\ \mathrm{d}\Omega 
   - \int_\Gamma \boldsymbol{\delta u} \cdot \boldsymbol{t}\ \mathrm{d}\Gamma 
   - \int_\Omega \boldsymbol{\delta u} \cdot \boldsymbol{b}\ \mathrm{d}\Omega = 0
```
*Implementation*
```julia
G = 80e3; K = 160e3; # Shear and bulk modulus (material parameters)
tn = 1.0, b=Vec((0.0, 0.0, -1.0)); # Normal traction and body force (external loading)
material = WeakForm((δu, ∇δu, u, ∇u, u_dot, ∇u_dot) -> (∇δu ⊡ (2*G*dev(symmetric(∇u)) + 3*K*vol(∇u)) )
nh = NeumannHandler(dh)
add!(nh, Neumann(:u, 2, getfaceset(dh.grid, "right"), (x,t,n)->tn*n))
add!(nh, BodyLoad(:c, 2, (x,t)->b))
```
"""
struct WeakForm{F<:Function}
    internal_contribution::F
end

function FerriteAssembly.element_residual!(re, state, ae, material::WeakForm, cv, dh_fh, Δt, buffer)
    ae_old = FerriteAssembly.get_aeold(buffer)
    for q_point in 1:getnquadpoints(cv)
        dΩ = getdetJdV(cv, q_point)
        u = function_value(cv, q_point, ae)
        ∇u = function_gradient(cv, q_point, ae)
        u_dot = (u - function_value(cv, q_point, ae_old))/Δt
        ∇u_dot = (∇u - function_gradient(cv, q_point, ae_old))/Δt
        for i in 1:getnbasefunctions(cv)
            δNᵢ = shape_value(cv, q_point, i)
            ∇δNᵢ = shape_gradient(cv, q_point, i)
            re[i] += material.internal_contribution(δNᵢ, ∇δNᵢ, u, ∇u, u_dot, ∇u_dot)*dΩ
        end
    end
end

