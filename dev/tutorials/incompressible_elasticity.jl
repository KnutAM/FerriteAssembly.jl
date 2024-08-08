using LinearAlgebra
using Ferrite, FerriteAssembly

struct LinearElasticity{T}
    G::T
    K::T
end
LinearElasticity(;E, ν) = LinearElasticity(E/(2(1 + ν)), E*ν/((1+ν)*(1-2ν)))

function FerriteAssembly.element_routine!(Ke, fe, state, ae, mp::LinearElasticity, cv, buffer)
    cvu = cv.u
    cvp = cv.p
    # Extract the dof range for each field
    dru = dof_range(buffer, :u)
    drp = dof_range(buffer, :p)

    for q_point in 1:getnquadpoints(cvu)
        dΩ = getdetJdV(cvu, q_point)

        for (iᵤ, Iᵤ) in pairs(dru)
            ∇δNu_dev = dev(shape_symmetric_gradient(cvu, q_point, iᵤ))
            for (jᵤ, Jᵤ) in pairs(dru)
                ∇Nu_dev = dev(shape_symmetric_gradient(cvu, q_point, jᵤ))
                Ke[Iᵤ,Jᵤ] += 2 * mp.G * ∇δNu_dev ⊡ ∇Nu_dev * dΩ
            end
            div_δNu = shape_divergence(cvu, q_point, iᵤ)
            for (jₚ, Jₚ) in pairs(drp)
                Np = shape_value(cvp, q_point, jₚ)
                Ke[Iᵤ,Jₚ] -= Np * div_δNu * dΩ
            end
        end
        for (iₚ, Iₚ) in pairs(drp)
            δNp = shape_value(cvp, q_point, iₚ)
            for (jᵤ, Jᵤ) in pairs(dru)
                div_Nu = shape_divergence(cvu, q_point, jᵤ)
                Ke[Iₚ,Jᵤ] -= δNp * div_Nu * dΩ
            end
            for (jₚ, Jₚ) in pairs(drp)
                Np = shape_value(cvp, q_point, jₚ)
                Ke[Iₚ,Jₚ] -= 1/mp.K * δNp * Np * dΩ
            end
        end
    end
end

function create_cook_grid(;nx, ny)
    corners = [Vec{2}((0.0,   0.0)),
               Vec{2}((48.0, 44.0)),
               Vec{2}((48.0, 60.0)),
               Vec{2}((0.0,  44.0))]
    return generate_grid(Triangle, (nx, ny), corners);
end;

function solve(;ν, ipu, ipp)
    # Material behavior
    mp = LinearElasticity(;E=1.0, ν=ν)

    # Grid and dofhandler
    grid = create_cook_grid(;nx=50, ny=50)
    dh = DofHandler(grid)
    add!(dh, :u, ipu) # displacement
    add!(dh, :p, ipp) # pressure
    close!(dh)

    # Boundary conditions
    ch = ConstraintHandler(dh)
    add!(ch, Dirichlet(:u, getfacetset(dh.grid, "left"), (x,t) -> zero(Vec{2})))
    close!(ch)
    update!(ch, 0.0)

    lh = LoadHandler(dh)
    add!(lh, Neumann(:u, 3, getfacetset(dh.grid, "right"), Returns(Vec{2}((0.0, 1/16)))))

    # Cellvalues
    qr = QuadratureRule{RefTriangle}(3)
    ip_geo = Lagrange{RefTriangle,1}()
    cvu = CellValues(qr, ipu, ip_geo)
    cvp = CellValues(qr, ipp, ip_geo)

    # Setup assembly
    cv = (u=cvu, p=cvp) # Will be replaced by MultiCellValues in the future
    domainbuffer = setup_domainbuffer(DomainSpec(dh, mp, cv))

    # Allocate and do the assembly
    K = allocate_matrix(dh)
    f = zeros(ndofs(dh))
    assembler = start_assemble(K, f)
    work!(assembler, domainbuffer) # Assemble the stiffness matrix
    apply!(f, lh, 0.0)             # Apply Neumann BC
    apply!(K, f, ch)              # Apply Dirichlet BC

    # Solve the equation system
    u = Symmetric(K) \ f;          # Solve the equation system

    # Export the results
    filename = "cook_" * (isa(ipu, Lagrange{2,RefTetrahedron,1}) ? "linear" : "quadratic") *
                         "_linear"
    VTKGridFile(filename, dh) do vtk
        write_solution(vtk, dh, u)
    end
    return u
end

linear_ip    = Lagrange{RefTriangle,1}()
quadratic_ip = Lagrange{RefTriangle,2}()

u1 = solve(;ν=0.4999999, ipu=linear_ip^2, ipp=linear_ip)
u2 = solve(;ν=0.4999999, ipu=quadratic_ip^2, ipp=linear_ip);

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
