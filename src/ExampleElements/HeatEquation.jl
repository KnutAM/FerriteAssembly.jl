@doc raw"""
    StationaryFourier(k)

For solving stationary linear heat conduction (which uses Fourier's law) with conductivity `k`, 
such that the heat flux is ``\boldsymbol{q}=-k \nabla T``, where ``T`` is the temperature field. 

The strong form is,
```math
 \nabla \cdot \boldsymbol{q} = h  \quad \textbf{x} \in \Omega,
```
and the corresponding weak form is 
```math
   \int_\Omega [\nabla \delta T] \cdot \boldsymbol{q} \mathrm{d}\Omega 
   = \int_\Gamma \delta T\ q_\mathrm{n} \mathrm{d}\Gamma 
   - \int_\Omega \delta T\ h \mathrm{d}\Omega
```
where, on the right hand side, ``q_\mathrm{n}`` is a heat flux normal to the boundary ``\Gamma``,
and ``h`` is a volumetric heat supply. 
These contributions are not included in the element, and should be added with `FerriteNeumann.jl`
"""
struct StationaryFourier{T}
    k::T # Thermal conductivity
end

function FerriteAssembly.element_routine!(Ke, re, state, ae, m::StationaryFourier, cv, args...)
    n_basefuncs = getnbasefunctions(cv)
    # Loop over quadrature points
    for q_point in 1:getnquadpoints(cv)
        dΩ = getdetJdV(cv, q_point)
        q = -m.k*function_gradient(cv, q_point, ae)
        for i in 1:n_basefuncs
            ∇δN = shape_gradient(cv, q_point, i)
            re[i] += (∇δN ⋅ q) * dΩ
            for j in 1:n_basefuncs
                ∇N = shape_gradient(cv, q_point, j)
                Ke[i, j] += -m.k*(∇δN ⋅ ∇N) * dΩ
            end
        end
    end
end

@doc raw"""
    TransientFourier(k, c)

For solving the transient linear heat equation (which uses Fourier's law) with conductivity `k`, 
such that the heat flux is ``\boldsymbol{q}=-k \nabla T``, where ``T`` is the temperature field. 

The strong form is,
```math
    c\ \dot{T} + \nabla \cdot \boldsymbol{q} = h  \quad \textbf{x} \in \Omega,
```
and the corresponding time-discretized weak form is 
```math
    \int_\Omega \delta T\ c\ \frac{T - {}^\mathrm{n}T}{\Delta t}\ \mathrm{d}\Omega
    + \int_\Omega [\nabla \delta T] \cdot \boldsymbol{q}\ \mathrm{d}\Omega 
   = \int_\Gamma \delta T\ q_\mathrm{n}\ \mathrm{d}\Gamma 
   - \int_\Omega \delta T\ h\ \mathrm{d}\Omega
```
where ``{}^\mathrm{n}T`` is the old temperature (in the previous timestep) and ``\Delta t`` is the timestep. 
On the right hand side, ``q_\mathrm{n}`` is a heat flux normal to the boundary ``\Gamma``,
and ``h`` is a volumetric heat supply. 
These external contributions on the right hand side are not included in the element, 
and should be added with `FerriteNeumann.jl`
"""
struct TransientFourier{T}
    k::T # Thermal conductivity
    c::T # Heat capacity
end

function FerriteAssembly.element_routine!(Ke, re, state, ae, m::TransientFourier, cv, dh_fh, Δt, buffer)
    n_basefuncs = getnbasefunctions(cv)
    ae_old = FerriteAssembly.get_aeold(buffer)
    for q_point in 1:getnquadpoints(cv)
        # Calculate values for the current quadrature point
        dΩ = getdetJdV(cv, q_point)
        ∇T = function_gradient(cv, q_point, ae)
        T = function_value(cv, q_point, ae) 
        Told = function_value(cv, q_point, ae_old)
        Tdot = (T-Told)/Δt
        q = -m.k*∇T
        # Assemble values into residual and stiffness arrays
        for i in 1:n_basefuncs
            δNi = shape_value(cv, q_point, i)
            ∇δNi = shape_gradient(cv, q_point, i)
            re[i] += (δNi*m.c*Tdot + ∇δNi ⋅ q)*dΩ
            for j in 1:n_basefuncs
                ∇Nj = shape_gradient(cv, q_point, j)
                Nj = shape_value(cv, q_point, j)
                Ke[i, j] += (m.c*δNi*Nj/Δt - m.k*(∇δNi ⋅ ∇Nj)) * dΩ
            end
        end
    end
end
