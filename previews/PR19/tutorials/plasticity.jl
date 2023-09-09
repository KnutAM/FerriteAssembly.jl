using Printf, LinearAlgebra
using Tensors, MaterialModelsBase, Ferrite, FerriteAssembly

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

struct J2PlasticityState{T} <: AbstractMaterialState
    ϵp::SymmetricTensor{2,3,T,6} # plastic strain
    κ::T                         # hardening stress
end;
function MaterialModelsBase.initial_material_state(::J2PlasticityNew)
    return J2PlasticityState(zero(SymmetricTensor{2,3}), 0.0)
end;

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

grid = generate_grid(Tetrahedron, (20,2,4), zero(Vec{3}), Vec((10.0,1.0,1.0)))

ip = Lagrange{3, RefTetrahedron, 1}()
dh = DofHandler(grid); add!(dh, :u, 3, ip); close!(dh)

ch = ConstraintHandler(dh)
add!(ch, Dirichlet(:u, getfaceset(grid, "left"), Returns(zero(Vec{3}))))
close!(ch)

cellvalues = CellVectorValues(QuadratureRule{3,RefTetrahedron}(2), ip)

buffer = setup_domainbuffer(DomainSpec(dh, material, cellvalues));

f = zeros(ndofs(dh)) # Pre-allocate external force vector to only apply external load once per time step.
traction_function(t) = 1e7*t
lh = LoadHandler(dh)
qr_order = 3
add!(lh, Neumann(:u, qr_order, getfaceset(grid, "right"), (x,t,n)->Vec((0.0, 0.0, traction_function(t)))))

time_history = collect(range(0.5, 1.0, 10))

function solve_problem(buffer, time_history, dh, ch, lh)
    K = create_sparsity_pattern(dh)
    r = zeros(ndofs(dh))
    a = zeros(ndofs(dh))
    t_old = 0.0
    u_max = [0.0]   # To save for plotting
    for t in time_history
        newton_itr = -1
        print("\n t = $t:\n")
        # Set the time increment passed to the material model
        # Not actually needed for this material, but in general it should be set
        set_time_increment!(buffer, t - t_old)

        # Update and apply Dirichlet boundary conditions
        update!(ch, t)   # Update the constraint handler and
        apply!(a, ch)    # set the prescribed values in the solution vector

        # Update and apply Neumann boundary conditions
        fill!(f, 0)      # Reset the external load
        apply!(f, lh, t) # Apply the Neumann boundary conditions
        while true; newton_itr += 1
            newton_itr > 8 && error("Reached maximum Newton iterations, aborting")
            # Assemble the contributions
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

        # Update the old states with the converged values for next timestep
        update_states!(buffer)  # I.e states_old = states
        t_old = t   # Update the old time (to calculate Δt)
        push!(u_max, maximum(abs, a)) # Save the maximum displacement in current timestep
    end
    return a, u_max
end
a, u_max = solve_problem(buffer, time_history, dh, ch, lh)

traction = vcat(0.0, traction_function.(time_history))

struct PlasticIntegrator{T}
    σvm::Vector{T}
    κ::Vector{T}
end
PlasticIntegrator(ncells::Int) = PlasticIntegrator(fill(NaN, ncells), fill(NaN, ncells))
intval = PlasticIntegrator(getncells(grid))
integrator = Integrator(intval)
function FerriteAssembly.integrate_cell!(intval::PlasticIntegrator, states, ae, m::J2PlasticityNew, cv, buffer)
    cellnr = cellid(buffer)
    V = σvm = κ = 0.0
    for q_point in 1:getnquadpoints(cv)
        state = states[q_point]
        ϵ = function_symmetric_gradient(cv, q_point, ae)
        σ = m.D⊡(ϵ - state.ϵp)
        dΩ = getdetJdV(cv, q_point)
        V += dΩ
        σvm += vonmises(σ)*dΩ
        κ += state.κ*dΩ
    end
    intval.σvm[cellnr] = σvm/V
    intval.κ[cellnr] = κ/V
end

work!(integrator, buffer; a=a)

vtk_grid("plasticity", dh) do vtkfile
    vtk_point_data(vtkfile, dh, a) # displacement field
    vtk_cell_data(vtkfile,  intval.σvm, "von Mises [Pa]")
    vtk_cell_data(vtkfile, intval.κ, "Drag stress [Pa]")
end

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
