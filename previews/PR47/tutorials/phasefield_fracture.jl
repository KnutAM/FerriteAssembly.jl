#=
# Phase-field fracture

TODO: Describe theory, use micromorphic variable from Ritu...
=#
using Ferrite, FerriteAssembly, FerriteMeshParser
using AppleAccelerate

## Implementation of physics

struct PhaseFieldFracture{C, T}
    G::T    # Elastic shear modulus
    K::T    # Elastic bulk modulus
    Gc::T   # Fracture energy
    l::T    # Length parameter
    α::T    # Micromorphic penalty factor
end
function PhaseFieldFracture(;E, ν, Gc, l, β)
    G = E / (2 * (1 + ν))
    K = E / (3 * (1 - 2ν))
    α = β * Gc / l
    T = promote_type(typeof(G), typeof(α))
    return PhaseFieldFracture{Nothing, T}(G, K, Gc, l, α)
end
PhaseFieldFracture{C}(m::PhaseFieldFracture{<:Any, T}) where {C, T} = PhaseFieldFracture{C, T}(m.G, m.K, m.Gc, m.l, m.α)

# Elastic (displacement) part (:u)
function FerriteAssembly.element_residual!(re, state, ae, m::PhaseFieldFracture{:u}, cv::CellValues, buffer)
    # Fixed parameters
    gϕ_min = 1e-10
    ## Values from phase-field problem
    cb_d = FerriteAssembly.get_coupled_buffer(buffer, :d)
    phasefields = FerriteAssembly.get_state(cb_d)

    for q_point in 1:getnquadpoints(cv)
        dΩ = getdetJdV(cv, q_point)
        ϵ = function_symmetric_gradient(cv, q_point, ae)
        ϕ = phasefields[q_point]
        gϕ = (1 - ϕ)^2 # AT1 and AT2
        gϕ_reg = (1 - gϕ_min)*gϕ + gϕ_min
        σ = (gϕ_reg * 2 * m.G) * dev(ϵ) + (gϕ_reg * m.K) * vol(ϵ)
        for i in 1:getnbasefunctions(cv)
            ∇δNu = shape_symmetric_gradient(cv, q_point, i)
            re[i] += (∇δNu ⊡ σ) * dΩ
        end
    end
end

# Phase-Field (damage) part (:d)
function FerriteAssembly.element_residual!(re, state, ae, m::PhaseFieldFracture{:d}, cv::CellValues, buffer)
    ## Values from elasticity problem
    cb_u = FerriteAssembly.get_coupled_buffer(buffer, :u)
    ae_u = FerriteAssembly.get_ae(cb_u)
    cv_u = FerriteAssembly.get_values(cb_u)
    
    ## Old phasefield values
    state_old = FerriteAssembly.get_old_state(buffer)
    for q_point in 1:getnquadpoints(cv)
        dΩ = getdetJdV(cv, q_point)
        d = function_value(cv, q_point, ae)
        ∇d = function_gradient(cv, q_point, ae)
        ⁿϕ = state_old[q_point]

        ## Elastic energy (no split)
        ϵ = function_symmetric_gradient(cv_u, q_point, ae_u)
        Ψ = 0.5 * (m.K - 2 * m.G / 3) * tr(ϵ)^2 + m.G * (ϵ ⊡ ϵ)

        ## AT1 model
        cw = 8/3
        ϕ = min(max((2 * Ψ + m.α * d - 3 * m.Gc/(8 * m.l))/(2 * Ψ + m.α), ⁿϕ), 1)
        ## AT2 model
        cw = 2.0
        ϕ = min(max((2 * Ψ + m.α * d)/(2 * Ψ + m.α + m.Gc / m.l), ⁿϕ), 1)
        
        for i in 1:getnbasefunctions(cv)
            ∇δNd = shape_gradient(cv, q_point, i)
            δNd = shape_value(cv, q_point, i)
            re[i] += ((2 * m.Gc * m.l / cw) * (∇δNd ⋅ ∇d) - m.α * (ϕ - d) * δNd) * dΩ
        end
        state[q_point] = FerriteAssembly.remove_dual(ϕ)
    end
end

function FerriteAssembly.create_cell_state(::PhaseFieldFracture{:d}, cv::CellValues, x, ae, args...)
    return [function_value(cv, i, ae) for i in 1:getnquadpoints(cv)]
end

## Simulation setup
# We start by loading the grid
grid = get_ferrite_grid(joinpath(@__DIR__, "sent_fine.inp"));
mbase = PhaseFieldFracture(;E = 210e3, ν = 0.3, Gc = 2.7, l = 1 * 1.5e-2, β = 100.0)

# Before defining the quadrature rules that are the same for both parts
qr_tri = QuadratureRule{RefTriangle}(2);
qr_quad = QuadratureRule{RefQuadrilateral}(2);

### Generic setup
function setup(m, grid, fieldname; qr_tri, qr_quad, ip_tri, ip_quad)
    dh = DofHandler(grid)
    
    sdh_tri = SubDofHandler(dh, getcellset(grid, "CPS3"))
    add!(sdh_tri, fieldname, ip_tri)
    cv_tri = CellValues(qr_tri, ip_tri)

    sdh_quad = SubDofHandler(dh, getcellset(grid, "CPS4R"))
    add!(sdh_quad, fieldname, ip_quad)
    cv_quad = CellValues(qr_quad, ip_quad)

    close!(dh)

    domains = Dict(
        "tri"  => DomainSpec(sdh_tri, m, cv_tri),
        "quad" => DomainSpec(sdh_quad, m, cv_quad))
    
    db = setup_domainbuffers(domains; threading = true, autodiffbuffer = true, a = zeros(ndofs(dh)))
    K = allocate_matrix(dh)
    r = zeros(ndofs(dh))
    return db, K, r, ndofs(dh)
end

db_u_uc, Ku, ru, ndofs_u = setup(PhaseFieldFracture{:u}(mbase), grid, :u; 
    qr_tri, qr_quad, 
    ip_tri = Lagrange{RefTriangle, 1}()^2, 
    ip_quad = Lagrange{RefQuadrilateral, 1}()^2
    )

db_d_uc, Kd, rd, ndofs_d = setup(PhaseFieldFracture{:d}(mbase), grid, :d; 
    qr_tri, qr_quad,
    ip_tri = Lagrange{RefTriangle, 1}(), 
    ip_quad = Lagrange{RefQuadrilateral, 1}()
    )

sim_u = Simulation(couple_buffers(db_u_uc; d = db_d_uc), zeros(ndofs_u), zeros(ndofs_u))
sim_d = Simulation(couple_buffers(db_d_uc; u = db_u_uc), zeros(ndofs_d), zeros(ndofs_d))

# Setup loading and boundary conditions
function load_function(t, rates = [45.0 => 1e-4, Inf => 1e-6])
    v0 = 0.0
    for i in eachindex(rates)
        t0 = i == 1 ? zero(t) : rates[i-1][1]
        if t ≤ rates[i][1]
            return v0 + rates[i][2] * (t - t0)
        else
            v0 += rates[i][2] * (rates[i][1] - t0)
        end
    end
end

ch_u = ConstraintHandler(FerriteAssembly.get_dofhandler(sim_u))
add!(ch_u, Dirichlet(:u, getfacetset(grid, "bottom"), Returns(zero(Vec{2}))))
add!(ch_u, Dirichlet(:u, getfacetset(grid, "top"), (x, t) -> load_function(t), [2]))
close!(ch_u)

# Write function to solve one simulation part, given the other as input. 
function solve_single_part(sim, coupled, K, r, ch; firsttol = 1e-4, tol = 1e-6, maxiter = 100)
    for i in 1:maxiter
        assembler = start_assemble(K, r)
        work!(assembler, sim, coupled)
        ch === nothing || apply_zero!(K, r, ch)
        res = norm(r)
        if i == 1 && res < firsttol
            return true  # no modification required, already converged
        elseif res < tol
            return false # current part was modified to converge
        elseif i ≥ maxiter
            error("single part iterations didn't converge")
        end
        sim.a .-= K \ r # Update unknowns
    end
end

function solve(sim_u, sim_d, Ku, ru, Kd, rd, ch_u, grid)
    time_vector = collect(1:140)
    for (n, t) in enumerate(time_vector)
        update!(ch_u, t)
        apply!(sim_u.a, ch_u)
        num = 0
        max_staggered = 2500
        for iter in 1:max_staggered
            num = iter
            u_converged = solve_single_part(sim_u, CoupledSimulations(d = sim_d), Ku, ru, ch_u)
            u_converged && break # Displacement was converged without updating
            d_converged = solve_single_part(sim_d, CoupledSimulations(u = sim_u), Kd, rd, nothing)
            d_converged && break # Damage was converged without updating
            iter ≥ max_staggered && error("Did not converge in staggered iterations")
        end
        println(n, ": ", num)
        update_states!(sim_d) # Only d has state variables
        copyto!(sim_d.aold, sim_d.a)
        copyto!(sim_u.aold, sim_u.a)
        VTKGridFile("fracture-$n", grid) do vtk
            write_solution(vtk, FerriteAssembly.get_dofhandler(sim_d), sim_d.a)
            write_solution(vtk, FerriteAssembly.get_dofhandler(sim_u), sim_u.a)
        end
    end
end
