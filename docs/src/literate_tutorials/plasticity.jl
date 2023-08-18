# # Plasticity
# This tutorial shows how to solve a plasticity problem by using the interface to  
# [`MaterialModelsBase.jl`](https://github.com/KnutAM/MaterialModelsBase.jl)
# included in `FerriteAssembly.jl`. Specifically, FerriteAssembly comes with the 
# `element_routine!` defined for any `MaterialModelsBase.AbstractMaterial`.
# 
# We adopt the material model as well as the simulation case from `Ferrite.jl`'s
# [plasticity example](https://ferrite-fem.github.io/Ferrite.jl/stable/examples/plasticity/)

using Printf, LinearAlgebra
using Tensors, MaterialModelsBase, Ferrite, FerriteAssembly

# ## Material modeling 
# We start by defining the J2PlasticityNew material model parameter struct, and an instance 
# thereof
struct J2PlasticityNew{T} <: AbstractMaterial
    G::T  # Shear modulus
    K::T  # Bulk modulus
    σ0::T # Initial yield limit
    H::T  # Hardening modulus
    D::SymmetricTensor{4, 3, T, 36} # Elastic stiffness tensor
end
function J2PlasticityNew(E, ν, σ₀, H)
    G = E / 2(1 + ν)
    K = E / 3(1 - 2ν)
    I2 = one(SymmetricTensor{2,3})
    I4 = one(SymmetricTensor{4,3})
    IxI = I2 ⊗ I2 
    D = 2G*(I4 - IxI/3) + K*IxI 
    return J2PlasticityNew(G, K, σ₀, H, D)
end
material = J2PlasticityNew(200.0e9, 0.3, 200.0e6, 10.0e9);

# Followed by the state variable struct along with its initial condition.
struct J2PlasticityState{T} <: AbstractMaterialState
    ϵp::SymmetricTensor{2,3,T,6} # plastic strain
    κ::T                         # hardening stress
end;
function MaterialModelsBase.initial_material_state(::J2PlasticityNew)
    return J2PlasticityState(zero(SymmetricTensor{2,3}), 0.0)
end;

# And finally the actual `material_response` function, with a few helper functions
function vonmises(σ)
    σdev = dev(σ)
    return sqrt((3/2)*(σdev⊡σdev))
end
function calculate_plastic_stress(ϵ, m, state) # When plastic
    σ_trial = m.D ⊡ (ϵ - state.ϵp)
    σeff_trial = vonmises(σ_trial)
    σdev_trial = dev(σ_trial)
    Φ_trial = σeff_trial - (m.σ0 + state.κ)
    μ = Φ_trial / (m.H + 3m.G)       # plastic multiplier
    σdev = (1 - 3m.G * μ / σeff_trial) * σdev_trial       
    return σdev + vol(σ_trial)   # updated stress
end
function MaterialModelsBase.material_response(
    material::J2PlasticityNew, ϵ::SymmetricTensor{2,3}, state::J2PlasticityState, 
    Δt, cache=get_cache(material), args...; kwargs...)

    σ_trial = material.D ⊡ (ϵ - state.ϵp) # trial-stress
    Φ_trial = vonmises(σ_trial) - (material.σ0 + state.κ)
    if Φ_trial < 0.0 # elastic loading
        return σ_trial, material.D, J2PlasticityState(state.ϵp, state.κ)
    else # plastic loading
        dσdϵ, σ = gradient(ϵ_ -> calculate_plastic_stress(ϵ_, material, state), ϵ, :all)

        μ =  Φ_trial / (material.H + 3*material.G)   # plastic multiplier
        κ = state.κ + μ*material.H
        σeff = material.σ0 + κ
        
        return σ, dσdϵ, J2PlasticityState(state.ϵp + (μ*3/(2*σeff))*dev(σ), κ)
    end
end;

# ## Standard `Ferrite.jl` setup
# With all required functions defined, we can now setup and assemble the finite element problem 
# using only Ferrite functionality
grid = generate_grid(Tetrahedron, (20,2,4), zero(Vec{3}), Vec((10.0,1.0,1.0)))

ip = Lagrange{3, RefTetrahedron, 1}()
dh = DofHandler(grid); add!(dh, :u, 3, ip); close!(dh)

ch = ConstraintHandler(dh) 
add!(ch, Dirichlet(:u, getfaceset(grid, "left"), Returns(zero(Vec{3}))))
close!(ch)

cellvalues = CellVectorValues(QuadratureRule{3,RefTetrahedron}(2), ip)

# ## Setting up the assembly
# Using the `setup_domainbuffer` function, 
buffer = setup_domainbuffer(DomainSpec(dh, material, cellvalues));
# we setup the `buffer`, old state variables, and new state variables. 
# The state variables are created via the [`create_cell_state`](@ref FerriteAssembly.create_cell_state) 
# function that is already defined for `MaterialModelsBase.AbstractMaterial`, since we overloaded 
# `MaterialModelsBase.initial_material_state` above. 
# 
# So far, we haven't included any loading, and following the original example, we would add loading 
# as a Neumann boundary condition on the right side. For this, we'll use FerriteAssembly's convenience 
# LoadHandler
f = zeros(ndofs(dh)) # Pre-allocate external force vector to only apply external load once per time step. 
traction_function(t) = 1e7*t
lh = LoadHandler(dh)
qr_order = 3
add!(lh, Neumann(:u, qr_order, getfaceset(grid, "right"), (x,t,n)->Vec((0.0, 0.0, traction_function(t)))))

# ## Solving the problem
# We first define how to step in time, 
# and following the original example this would be
time_history = collect(range(0.5, 1.0, 10))

# We can now solve the problem using Newton iterations, for each time step.
function solve_problem(buffer, time_history, dh, ch, lh)
    K = create_sparsity_pattern(dh)
    r = zeros(ndofs(dh))
    a = zeros(ndofs(dh))
    t_old = 0.0
    u_max = [0.0]   # To save for plotting
    for t in time_history
        newton_itr = -1
        print("\n t = $t:\n")
        ## Set the time increment passed to the material model
        ## Not actually needed for this material, but in general it should be set
        set_time_increment!(buffer, t - t_old)

        ## Update and apply Dirichlet boundary conditions
        update!(ch, t)   # Update the constraint handler and
        apply!(a, ch)    # set the prescribed values in the solution vector

        ## Update and apply Neumann boundary conditions
        fill!(f, 0)      # Reset the external load
        apply!(f, lh, t) # Apply the Neumann boundary conditions
        while true; newton_itr += 1
            newton_itr > 8 && error("Reached maximum Newton iterations, aborting")
            ## Assemble the contributions
            assembler = start_assemble(K, r) # K and r are zeroed by this call. 
            work!(assembler, buffer; a=a)
            r .-= f                 # Subtract external forces
            apply_zero!(K, r, ch)   # Apply Dirichlet boundary conditions
            
            norm_r = norm(r) # Constrained dofs in r where zeroed by apply_zero!
            print("Iteration: $newton_itr \tresidual: $(@sprintf("%.8f", norm_r))\n")
            norm_r < 1.0 && break # Tolerance 1 N

            Δa = -Symmetric(K) \ r  # Do one Newton update
            apply_zero!(Δa, ch)     # Ensure exact BC (See Ferrite's doc)
            a .+= Δa                # Update the displacements
        end

        ## Update the old states with the converged values for next timestep
        update_states!(buffer)  # I.e states_old = states
        t_old = t   # Update the old time (to calculate Δt)
        push!(u_max, maximum(abs, a)) # Save the maximum displacement in current timestep
    end
    return a, u_max
end
a, u_max = solve_problem(buffer, time_history, dh, ch, lh)

# ## Postprocessing 
# As in the Ferrite example, we would like to plot the traction versus the 
# maximum displacements, and export the final displacement field, von Mises
# stress, and hardening stress. See Ferrite's example for the plotting. 
traction = vcat(0.0, traction_function.(time_history))

# Instead of using the average of the integration point values, as in Ferrite,
# we'll perform the volume averaging (which for 2nd order integration is equivalent),
# by using the Integrator 

## test the result                       #src
using Test                               #src
@test norm(u_max[end]) ≈ 0.254452645     #src
