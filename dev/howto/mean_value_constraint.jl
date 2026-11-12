using Ferrite, FerriteAssembly
import FerriteAssembly.ExampleElements: StationaryFourier

struct MeanValueHeat{T}
    k::T # Conductivity
    f::T # Heat source
    ū::T # Mean temperature
end

function FerriteAssembly.element_routine!(Ke, re, state, ae, m::MeanValueHeat, cv, buffer)
    dofs_u = dof_range(buffer, :u)
    dofs_λ = dof_range(buffer, :λ)
    dof_λ = only(dofs_λ) # λ is a scalar and has one dof
    λ = algebraic_value(FerriteAssembly.get_algebraic_values(buffer, :λ), ae, dofs_λ)
    for q_point in 1:getnquadpoints(cv)
        dΩ = getdetJdV(cv, q_point)
        u = function_value(cv, q_point, ae, dofs_u)
        ∇u = function_gradient(cv, q_point, ae, dofs_u)
        re[dof_λ] += (u - m.ū) * dΩ
        for (i, I) in pairs(dofs_u)
            δN = shape_value(cv, q_point, i)
            ∇δN = shape_gradient(cv, q_point, i)
            re[I] += (m.k * (∇δN ⋅ ∇u) + δN * (λ - m.f)) * dΩ
            Ke[I, dof_λ] += δN * dΩ
            Ke[dof_λ, I] += δN * dΩ
            for (j, J) in pairs(dofs_u)
                ∇N = shape_gradient(cv, q_point, j)
                Ke[I, J] += m.k * (∇δN ⋅ ∇N) * dΩ
            end
        end
    end
end;

grid = generate_grid(Quadrilateral, (20, 20))
ip = Lagrange{RefQuadrilateral, 1}()
dh = DofHandler(grid)
add!(dh, :u, ip)
add!(dh, :λ, AlgebraicVariable())
close!(dh)
cv = CellValues(QuadratureRule{RefQuadrilateral}(2), ip);

material = MeanValueHeat(1.0, 0.5, 2.0)
buffer = setup_domainbuffer(DomainSpec(dh, material, cv; algebraic_per_cell = (:λ,)));

K = allocate_matrix(dh; algebraic_couplings = FerriteAssembly.algebraic_couplings(buffer))
r = zeros(ndofs(dh))
a = zeros(ndofs(dh));

work!(KeReAssembler(K, r), buffer; a = a)
a .-= K \ r;

algebraic_value(dh, a, :λ) ≈ material.f

integrator = SimpleIntegrator((u, ∇u, state) -> (1.0, u), (0.0, 0.0))
work!(integrator, buffer; a = a)
integrator.val[2] / integrator.val[1] ≈ material.ū

struct BoundaryMeanValue{T}
    ū::T
end

function FerriteAssembly.facet_routine!(Ke, re, ae, m::BoundaryMeanValue, fv, buffer)
    dofs_u = dof_range(buffer, :u)
    dofs_λ = dof_range(buffer, :λ)
    dof_λ = only(dofs_λ)
    λ = algebraic_value(FerriteAssembly.get_algebraic_values(buffer, :λ), ae, dofs_λ)
    for q_point in 1:getnquadpoints(fv)
        dΓ = getdetJdV(fv, q_point)
        u = function_value(fv, q_point, ae, dofs_u)
        re[dof_λ] += (u - m.ū) * dΓ
        for (i, I) in pairs(dofs_u)
            δN = shape_value(fv, q_point, i)
            re[I] += δN * λ * dΓ
            Ke[I, dof_λ] += δN * dΓ
            Ke[dof_λ, I] += δN * dΓ
        end
    end
end;

cellbuffer = setup_domainbuffer(DomainSpec(dh, StationaryFourier(material.k), cv))
lh = LoadHandler(dh)
add!(lh, BodyLoad(:u, 2, Returns(material.f)))

fv = FacetValues(FacetQuadratureRule{RefQuadrilateral}(2), ip)
boundary_material = BoundaryMeanValue(1.0)
facetbuffer = setup_domainbuffer(DomainSpec(dh, boundary_material, fv;
    set = getfacetset(grid, "right"), algebraic_per_cell = (:λ,)));

K = allocate_matrix(dh; algebraic_couplings = FerriteAssembly.algebraic_couplings(cellbuffer, facetbuffer))
fill!(a, 0)
work!(KeReAssembler(K, r), cellbuffer; a = a)
work!(KeReAssembler(K, r; fillzero = false), facetbuffer; a = a) # Don't reset K and r
fext = zeros(ndofs(dh))
apply!(fext, lh, 0.0)
a .-= K \ (r - fext);

integrator = SimpleIntegrator((u, ∇u, n) -> (1.0, u), (0.0, 0.0))
work!(integrator, facetbuffer; a = a)
integrator.val[2] / integrator.val[1] ≈ boundary_material.ū

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
