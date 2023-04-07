"""
    FourierMaterial(k)

For solving the linear heat equation (which uses Fourier's law) according to the strong form,
```math
 -\\nabla \\cdot (k \\nabla T) = f  \\quad \\textbf{x} \\in \\Omega,
```
where ``T`` is the temperature field, ``k`` the thermal conductivity, and ``f`` a volume source term.
``f`` is not included in the element assembly, and can be applied using `FerriteNeumann.jl` 
Only the stiffness contribution is assembled.
"""
struct FourierMaterial{T}
    k::T # Thermal conductivity
end

function FA.element_routine!(Ke, re, state, ae, m::FourierMaterial, cv, args...)
    n_basefuncs = getnbasefunctions(cv)
    # Loop over quadrature points
    for q_point in 1:getnquadpoints(cv)
        dΩ = getdetJdV(cv, q_point)
        for i in 1:n_basefuncs
            ∇δN = shape_gradient(cv, q_point, i)
            for j in 1:n_basefuncs
                ∇N = shape_gradient(cv, q_point, j)
                Ke[i, j] += m.k*(∇δN ⋅ ∇N) * dΩ
            end
        end
    end
end

"""
    TransientFourierMaterial(k, c)

For solving transient heat flow, assuming flux according to Fourier's law
with conductivity, `k`, and heat capacity, `c`, according to 
```math
  c\\frac{\\partial T}{\\partial t}-\\nabla \\cdot (k \\nabla T) = f  \\quad x \\in \\Omega,
```
where ``T`` is the temperature field and ``f`` is a volume source term. 
``f`` is not included in the element assembly, but can be included with `FerriteNeumann.jl`.
Both the stiffness matrix, ``K=\\mathrm{d}f_{int}/\\partial T`` 
and the residual vector, ``r=f_{int}-f_{ext}`` are assembled. 
"""
struct TransientFourierMaterial{T}
    k::T # Thermal conductivity
    c::T # Heat capacity
end

function FA.element_routine!(Ke, re, state, ae, m::TransientFourierMaterial, cv, dh_fh, Δt, buffer)
    n_basefuncs = getnbasefunctions(cv)
    ae_old = FA.get_aeold(buffer)
    for q_point in 1:getnquadpoints(cv)
        # Calculate values for the current quadrature point
        dΩ = getdetJdV(cv, q_point)
        ∇T = function_gradient(cv, q_point, ae)
        T = function_value(cv, q_point, ae) 
        Told = function_value(cv, q_point, ae_old)
        Tdot = (T-Told)/Δt
        # Assemble values into residual and stiffness arrays
        for i in 1:n_basefuncs
            δNi = shape_value(cv, q_point, i)
            ∇δNi = shape_gradient(cv, q_point, i)
            re[i] += (m.k*(∇δNi ⋅ ∇T) + δNi*m.c*Tdot)*dΩ
            for j in 1:n_basefuncs
                ∇Nj = shape_gradient(cv, q_point, j)
                Nj = shape_value(cv, q_point, j)
                Ke[i, j] += (m.k*(∇δNi ⋅ ∇Nj) + m.c*δNi*Nj ) * dΩ
            end
        end
    end
end
