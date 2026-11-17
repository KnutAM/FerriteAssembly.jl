using Ferrite, FerriteAssembly
import MaterialModelsBase as MMB
import MechanicalMaterialModels as MMM
import CairoMakie as CM

struct OutOfPlaneCellValues{CV <: CellValues} <: Ferrite.AbstractCellValues
    cv::CV                      # CellValues for the in-plane problem
    x::Vector{Vec{3, Float64}}  # Cell coordinates lifted to 3d
end
function OutOfPlaneCellValues(cv::CellValues)
    return OutOfPlaneCellValues(cv, zeros(Vec{3, Float64}, Ferrite.getngeobasefunctions(cv)))
end

# Lifting of 2d tensors to 3d
lift3d(v::Number) = v
lift3d(v::Vec{2, T}) where {T} = Vec{3, T}((v[1], v[2], zero(T)))
lift3d(t::Tensor{2, 2, T}) where {T} = Tensor{2, 3, T}((i, j) -> (i ≤ 2 && j ≤ 2) ? t[i, j] : zero(T))

# Gradient of the out-of-plane base function, for scalar and vector fields
out_of_plane_gradient(::Number) = Vec{3}((0.0, 0.0, 1.0))
out_of_plane_gradient(::Vec) = Tensor{2, 3}((i, j) -> (i == 3 && j == 3) ? 1.0 : 0.0)

# The required accessors
Ferrite.getnbasefunctions(v::OutOfPlaneCellValues) = getnbasefunctions(v.cv) + 1
Ferrite.getnquadpoints(v::OutOfPlaneCellValues) = getnquadpoints(v.cv)
Ferrite.getdetJdV(v::OutOfPlaneCellValues, q_point::Int) = getdetJdV(v.cv, q_point)
Ferrite.getngeobasefunctions(v::OutOfPlaneCellValues) = Ferrite.getngeobasefunctions(v.cv)
Ferrite.geometric_value(v::OutOfPlaneCellValues, q_point::Int, i::Int) = Ferrite.geometric_value(v.cv, q_point, i)

function Ferrite.reinit!(v::OutOfPlaneCellValues, cell::Union{Nothing, Ferrite.AbstractCell}, x::AbstractVector{<:Vec})
    v.x .= lift3d.(x)
    return reinit!(v.cv, cell, x)
end

function Ferrite.shape_value(v::OutOfPlaneCellValues, q_point::Int, i::Int)
    N1 = lift3d(shape_value(v.cv, q_point, 1)) # Ensure the same type for all base functions
    return i ≤ getnbasefunctions(v.cv) ? lift3d(shape_value(v.cv, q_point, i)) : zero(N1)
end

function Ferrite.shape_gradient(v::OutOfPlaneCellValues, q_point::Int, i::Int)
    i ≤ getnbasefunctions(v.cv) && return lift3d(shape_gradient(v.cv, q_point, i))
    return out_of_plane_gradient(shape_value(v.cv, q_point, 1))
end

FerriteAssembly.routine_coordinates(v::OutOfPlaneCellValues, ::AbstractVector{<:Vec}) = v.x;

qr = QuadratureRule{RefQuadrilateral}(2)
ip = Lagrange{RefQuadrilateral, 1}()^2
cv_oop = OutOfPlaneCellValues(CellValues(qr, ip))
let grid = generate_grid(Quadrilateral, (1, 1))
    FerriteAssembly.check_values_interface(cv_oop, getcells(grid, 1), getcoordinates(grid, 1))
end

L, H = (10.0, 2.0)
grid = generate_grid(Quadrilateral, (40, 8), Vec((0.0, 0.0)), Vec((L, H)))
transform_coordinates!(grid, x -> Vec((x[1], x[2] * (1 - 0.5 * x[1] / L))))

dh = DofHandler(grid)
add!(dh, :u, ip)
add!(dh, :ϵ33, AlgebraicVariable())
close!(dh);

function create_constrainthandler(dh, Δu; plane_strain = false)
    ch = ConstraintHandler(dh)
    add!(ch, Dirichlet(:u, getfacetset(dh.grid, "left"), Returns(0.0), 1))
    add!(ch, Dirichlet(:u, getfacetset(dh.grid, "bottom"), Returns(0.0), 2))
    add!(ch, Dirichlet(:u, getfacetset(dh.grid, "right"), (x, t) -> Δu * t, 1))
    if plane_strain
        add!(ch, AffineConstraint(only(algebraic_dofs(dh, :ϵ33)), Pair{Int, Float64}[], 0.0))
    end
    close!(ch)
    return ch
end;

function right_edge_dofs(dh)
    ch = ConstraintHandler(dh)
    add!(ch, Dirichlet(:u, getfacetset(dh.grid, "right"), Returns(0.0), 1))
    close!(ch)
    return ch.prescribed_dofs
end;

material = MMM.Plastic(;
    elastic = MMM.LinearElastic(; E = 210.0e3, ν = 0.3),
    yield = 300.0,
    isotropic = MMM.Voce(; Hiso = 10.0e3, κ∞ = 200.0),
    kinematic = MMM.ArmstrongFrederick(; Hkin = 20.0e3, β∞ = 100.0),
);

domain_gps = DomainSpec(dh, material, cv_oop;
    algebraic_variables = (FerriteAssembly.AttachAlgebraic(:ϵ33 => :u),));

function solve(domain, dh; ϵ33_dofs = Int[], plane_strain = false, Δu = 0.1, num_steps = 40, tolerance = 1.0e-8, maxiter = 20)
    buffer = setup_domainbuffer(domain)
    ch = create_constrainthandler(dh, Δu; plane_strain)
    right_dofs = right_edge_dofs(dh)
    K = allocate_matrix(dh, ch; algebraic_couplings = FerriteAssembly.algebraic_couplings(buffer))
    r = zeros(ndofs(dh))
    a = zeros(ndofs(dh))
    displacement = [0.0]
    force = [0.0]
    ϵ33 = [0.0]
    t_old = 0.0
    for t in range(0, 1, num_steps + 1)[2:end]
        update!(ch, t - t_old)
        reaction = 0.0
        for i in 1:maxiter
            assembler = start_assemble(K, r)
            work!(assembler, buffer; a)
            reaction = sum(r[right_dofs])
            if i == 1 # Apply the increments of the prescribed values
                r .*= -1
                apply!(K, r, ch)
                a .+= K \ r
            else
                apply_zero!(K, r, ch)
                norm(r) < tolerance * max(1.0, abs(reaction)) && break
                i == maxiter && error("Did not converge")
                a .-= K \ r
            end
        end
        update_states!(buffer)
        push!(displacement, Δu * t)
        push!(force, reaction)
        push!(ϵ33, isempty(ϵ33_dofs) ? NaN : a[only(ϵ33_dofs)])
        t_old = t
    end
    return displacement, force, ϵ33
end;

displacement, force_gps, ϵ33_gps = solve(domain_gps, dh; ϵ33_dofs = algebraic_dofs(dh, :ϵ33));

_, force_pe, _ = solve(domain_gps, dh; plane_strain = true);

dh_ps = DofHandler(grid)
add!(dh_ps, :u, ip)
close!(dh_ps)
domain_ps = DomainSpec(dh_ps, MMB.ReducedStressState(MMB.PlaneStress(), material), CellValues(qr, ip))
_, force_ps, _ = solve(domain_ps, dh_ps);

all(force_ps[2:end] .< force_gps[2:end] .< force_pe[2:end])

fig = CM.Figure()
ax_f = CM.Axis(fig[1, 1]; xlabel = "Displacement [mm]", ylabel = "Force per thickness [N/mm]")
CM.lines!(ax_f, displacement, force_pe; label = "Plane strain")
CM.lines!(ax_f, displacement, force_gps; label = "Generalized plane stress")
CM.lines!(ax_f, displacement, force_ps; label = "Plane stress", linestyle = :dash)
CM.axislegend(ax_f; position = :rb)
ax_e = CM.Axis(fig[2, 1]; xlabel = "Displacement [mm]", ylabel = "ϵ₃₃ [-]")
CM.lines!(ax_e, displacement, ϵ33_gps)
fig

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
