@doc raw"""
    WeakForm(f::Function)

Solve the problem with primary variable, ``u``, variation, ``\delta u``, and a weak form
```math
   \int_\Omega f(\delta u, \delta u \otimes\nabla, u, u \otimes\nabla, \dot{u}, \dot{u} \otimes\nabla) \mathrm{d}\Omega 
   - \int_\Gamma \delta u\ h(x,t,n)\ \mathrm{d}\Gamma - \int_\Omega \delta u\ b(x,t)\ \mathrm{d}\Omega = 0
```
where the function `f` is given to the weak form, and `h` and `b` are given via the [`LoadHandler`](@ref FerriteAssembly.LoadHandler).
Note the sign convention: `LoadHandler`'s `apply!` *adds* the value of the given function to
the residual/external-load vector, so to subtract a term as shown above, the corresponding
function must return the negative of that term (e.g. `-h` and `-b`), as demonstrated below.

!!! note "This element is intended for testing"
    It is not optimized for speed

# Examples 
## Transient heat flow
*Weak form*
```math
    \int_\Omega \delta u\ c\ \dot{u} + k\ [\nabla \delta u] \cdot [\nabla u]\ \mathrm{d}\Omega
   - \int_\Gamma \delta u\ q_\mathrm{n}\ \mathrm{d}\Gamma - \int_\Omega \delta u\ b\ \mathrm{d}\Omega = 0
```
where ``q_\mathrm{n}`` is the outward heat flux on the boundary.
*Implementation*
This implementation is equivalent to [`TransientFourier`](@ref FerriteAssembly.ExampleElements.TransientFourier),
but it is also possible to add the body load directly in the weak form (if desired).
```julia
c = 1.0; k = 1.0; # heat capacity and heat conductivity (material parameters)
qn = 1.0; b=1.0;  # Outward boundary heat flux and internal heat source (external loading)
material = WeakForm((δu, ∇δu, u, ∇u, u_dot, ∇u_dot) -> δu*c*u_dot + k*(∇δu ⋅ ∇u))
lh = LoadHandler(dh)
add!(lh, Neumann(:u, 2, getfacetset(dh.grid, "right"), (x,t,n)->-qn)) # rᵢ -= ∫ δuᵢ*qₙ dΓ
add!(lh, BodyLoad(:u, 1, (x,t)->-b))                                  # rᵢ -= ∫ δuᵢ*b dΩ
```

## Linear elasticity
This implementation is equivalent to [`ElasticPlaneStrain`](@ref FerriteAssembly.ExampleElements.ElasticPlaneStrain),
but it is also possible to add the body load directly in the weak form (if desired). 
*Weak form*
```math
   \int_\Omega [\boldsymbol{\delta u}\otimes\nabla]^\mathrm{sym} : \boldsymbol{\sigma}\ \mathrm{d}\Omega 
   - \int_\Gamma \boldsymbol{\delta u} \cdot \boldsymbol{t}\ \mathrm{d}\Gamma 
   - \int_\Omega \boldsymbol{\delta u} \cdot \boldsymbol{b}\ \mathrm{d}\Omega = 0
```
*Implementation*
```julia
G = 80e3; K = 160e3; # Shear and bulk modulus (material parameters)
tn = 1.0; b = Vec((0.0, 0.0, -1.0)); # Outward normal traction and body force (external loading)
material = WeakForm((δu, ∇δu, u, ∇u, u_dot, ∇u_dot) -> (∇δu ⊡ (2*G*dev(symmetric(∇u)) + 3*K*vol(∇u))))
lh = LoadHandler(dh)
add!(lh, Neumann(:u, 2, getfacetset(dh.grid, "right"), (x,t,n)->-tn*n)) # rᵢ -= ∫ δuᵢ⋅(tₙn) dΓ
add!(lh, BodyLoad(:u, 2, (x,t)->-b))                                    # rᵢ -= ∫ δuᵢ⋅b dΩ
```
"""
struct WeakForm{F<:Function}
    internal_contribution::F
end

function FerriteAssembly.element_residual!(re, state, ae, material::WeakForm, cv, buffer)
    ae_old = FerriteAssembly.get_aeold(buffer)
    Δt = FerriteAssembly.get_time_increment(buffer)
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