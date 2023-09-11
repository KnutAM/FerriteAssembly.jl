using Ferrite, Tensors
using FerriteAssembly
import CairoMakie as CM

Base.@kwdef struct ZenerMaterial{T}
    K ::T=5.0   # Bulk modulus
    G1::T=1.0   # Shear modulus, parallel
    G2::T=50.   # Shear modulus, series
    η ::T=5.0   # Damping modulus
end;

function FerriteAssembly.create_cell_state(::ZenerMaterial, cv::CellVectorValues{dim}, args...) where dim
    return [zero(SymmetricTensor{2,dim}) for _ in 1:getnquadpoints(cv)]
end;

function FerriteAssembly.element_residual!(re, state, ae, m::ZenerMaterial, cv::CellValues, buffer)
    Δt = FerriteAssembly.get_time_increment(buffer)
    old_ϵvs = FerriteAssembly.get_old_state(buffer)
    for q_point in 1:getnquadpoints(cv)
        old_ϵv = old_ϵvs[q_point]
        dΩ = getdetJdV(cv, q_point)
        ϵ = function_symmetric_gradient(cv, q_point, ae)
        ϵdev = dev(ϵ)
        ϵv = (Δt*2*m.G2*ϵdev + m.η*old_ϵv)/(m.η + Δt*2*m.G2)
        σ = (m.G1+m.G2)*2*ϵdev - 2*m.G2*ϵv + m.K*vol(ϵ)
        for i in 1:getnbasefunctions(cv)
            δ∇N = shape_symmetric_gradient(cv, q_point, i)
            re[i] += (δ∇N⊡σ)*dΩ
        end
        # Note that to save the state by mutation, we need to extract the value from the dual
        # number. Consequently, we do this before assigning to the state vector. Note that
        # if the state was a scalar, we should use `ForwardDiff.value` instead.
        state[q_point] = Tensors._extract_value(ϵv)
    end
end;

grid = generate_grid(Quadrilateral, (2,2))
ip = Ferrite.default_interpolation(Quadrilateral)
dh = DofHandler(grid)
add!(dh, :u, 2, ip)
close!(dh)
qr = QuadratureRule{2,RefCube}(2)
cv = CellVectorValues(qr, ip, ip)
m = ZenerMaterial()
domain = DomainSpec(dh, m, cv)
buffer = setup_domainbuffer(domain);

ch = ConstraintHandler(dh)
add!(ch, Dirichlet(:u, getfaceset(grid, "left"), Returns(0.0), 1))
add!(ch, Dirichlet(:u, getfaceset(grid, "bottom"), Returns(0.0), 2))
close!(ch)
update!(ch, 0.0);

lh = LoadHandler(dh)
traction(t) = clamp(t, 0, 1)*Vec((1.0, 0.0))
add!(lh, Neumann(:u, 2, getfaceset(grid, "right"), (x, t, n) -> traction(t)));

function solve_nonlinear_timehistory(buffer, dh, ch, lh; time_history)
    maxiter = 10
    tolerance = 1e-6
    K = create_sparsity_pattern(dh)
    r = zeros(ndofs(dh))
    f = zeros(ndofs(dh))
    a = zeros(ndofs(dh))
    t_force = [0.0]
    u1_max = [0.0]
    told = 0.0
    for t in time_history
        # Update and apply the Dirichlet boundary conditions
        update!(ch, t)
        apply!(a, ch)
        # Update and apply the Neumann boundary conditions
        fill!(f, 0)
        apply!(f, lh, t)
        # Update the time increment (passed to `element_residual!`)
        set_time_increment!(buffer, t-told)
        for i in 1:maxiter
            # Assemble the system
            assembler = start_assemble(K, r)
            work!(assembler, buffer; a=a)
            r .-= f
            # Apply boundary conditions
            apply_zero!(K, r, ch)
            # Check convergence
            norm(r) < tolerance && break
            i == maxiter && error("Did not converge")
            # Solve the linear system and update the dof vector
            a .-= K\r
            apply_zero!(a, ch)
        end
        # If converged, update the old state variables to the current.
        update_states!(buffer)

        # Save values for postprocessing
        push!(t_force, norm(traction(t)))
        push!(u1_max, maximum(a))
        told = t
    end
    return u1_max, t_force
end;

time_history = collect(range(0,1,10)).^2
append!(time_history, 1 .+ collect(range(0,1,10)[2:end]).^2)

u1_max, t_force = solve_nonlinear_timehistory(buffer, dh, ch, lh; time_history=time_history[2:end]);

fig = CM.Figure()
ax_t = CM.Axis(fig[1,1]; xlabel="time", ylabel="traction")
ax_d = CM.Axis(fig[2,1]; xlabel="time", ylabel="displacement")
CM.lines!(ax_t, time_history, t_force)
CM.scatter!(ax_t, time_history, t_force)
CM.lines!(ax_d, time_history, u1_max)
CM.scatter!(ax_d, time_history, u1_max)
fig

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
