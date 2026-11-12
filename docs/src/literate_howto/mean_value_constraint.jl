# # Mean value constraints
# In this "how-to", we use `Ferrite`'s algebraic variables to prescribe the mean value of a field,
# both in the domain and on a part of the boundary. An algebraic variable is a global unknown that
# is not associated with the grid, here a Lagrange multiplier, $\lambda$, that enforces the constraint.
# The support for algebraic variables in `FerriteAssembly` is experimental, see
# [Algebraic variables](@ref) for further details.
#
# ## Mean value in the domain
# We consider the stationary heat equation with a heat source, $f$, where the heat flux is zero on the
# entire boundary. The temperature, $u$, is then only determined up to a constant. To obtain a unique
# solution, we prescribe the mean temperature, $\bar{u}$, giving the weak form
# ```math
# \begin{aligned}
# \int_\Omega \left[ k\, \nabla \delta u \cdot \nabla u - \delta u\, f + \delta u\, \lambda \right] \mathrm{d}\Omega &= 0 \\
# \delta \lambda \int_\Omega \left[ u - \bar{u} \right] \mathrm{d}\Omega &= 0
# \end{aligned}
# ```
# We define a material, `MeanValueHeat`, and the corresponding `element_routine!`. The local degrees
# of freedom consist of the cell's degrees of freedom followed by the degree of freedom for $\lambda$.
# `dof_range(buffer, name)` gives the local dof range, both for the field `:u` and the algebraic variable `:λ`.
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

# Note that `ae` contains the value of $\lambda$ after the values for $u$, so the dof range must
# be passed to `function_value` and `function_gradient`.
#
# We add the algebraic variable to the dofhandler, in addition to the temperature field.
grid = generate_grid(Quadrilateral, (20, 20))
ip = Lagrange{RefQuadrilateral, 1}()
dh = DofHandler(grid)
add!(dh, :u, ip)
add!(dh, :λ, AlgebraicVariable())
close!(dh)
cv = CellValues(QuadratureRule{RefQuadrilateral}(2), ip);

# When setting up the domain, we give the algebraic variables that are used in the domain
material = MeanValueHeat(1.0, 0.5, 2.0)
buffer = setup_domainbuffer(DomainSpec(dh, material, cv; algebraic_per_cell = (:λ,)));

# The sparsity pattern must include the entries for $\lambda$. The required information
# is created from the domain buffer and passed to `allocate_matrix`.
K = allocate_matrix(dh; algebraic_couplings = FerriteAssembly.algebraic_couplings(buffer))
r = zeros(ndofs(dh))
a = zeros(ndofs(dh));

# We can then assemble and solve the linear problem as usual
work!(KeReAssembler(K, r), buffer; a = a)
a .-= K \ r;

# As the heat flux is zero on the boundary, $\lambda$ must balance the heat source,
algebraic_value(dh, a, :λ) ≈ material.f

# and we can check that the mean temperature is as prescribed by integrating the solution
integrator = SimpleIntegrator((u, ∇u, state) -> (1.0, u), (0.0, 0.0))
work!(integrator, buffer; a = a)
integrator.val[2] / integrator.val[1] ≈ material.ū

# ## Mean value on the boundary
# Algebraic variables can also be used in facet domains. Here, we prescribe the mean
# temperature on the right boundary, $\Gamma$, instead of in the domain, i.e. we add the following
# contributions to the weak form,
# ```math
# \int_\Gamma \delta u\, \lambda\ \mathrm{d}\Gamma, \quad
# \delta \lambda \int_\Gamma \left[ u - \bar{u} \right] \mathrm{d}\Gamma
# ```
# which we implement in a `facet_routine!`,
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

# For the cells, we can now use the regular heat equation element, `StationaryFourier`, since the
# algebraic variable is only used in the facet domain. The heat source is added with a `LoadHandler`.
cellbuffer = setup_domainbuffer(DomainSpec(dh, StationaryFourier(material.k), cv))
lh = LoadHandler(dh)
add!(lh, BodyLoad(:u, 2, Returns(material.f)))

fv = FacetValues(FacetQuadratureRule{RefQuadrilateral}(2), ip)
boundary_material = BoundaryMeanValue(1.0)
facetbuffer = setup_domainbuffer(DomainSpec(dh, boundary_material, fv;
    set = getfacetset(grid, "right"), algebraic_per_cell = (:λ,)));

# The sparsity pattern is now created based on the facet domain. If algebraic variables are used
# in several domain buffers, all can be passed to `FerriteAssembly.algebraic_couplings`.
K = allocate_matrix(dh; algebraic_couplings = FerriteAssembly.algebraic_couplings(cellbuffer, facetbuffer))
fill!(a, 0)
work!(KeReAssembler(K, r), cellbuffer; a = a)
work!(KeReAssembler(K, r; fillzero = false), facetbuffer; a = a) # Don't reset K and r
fext = zeros(ndofs(dh))
apply!(fext, lh, 0.0)
a .-= K \ (r - fext);

# The mean temperature on the boundary is then as prescribed
integrator = SimpleIntegrator((u, ∇u, n) -> (1.0, u), (0.0, 0.0))
work!(integrator, facetbuffer; a = a)
integrator.val[2] / integrator.val[1] ≈ boundary_material.ū

# ## Multithreading
# The degree of freedom for $\lambda$ is shared between all cells (or facets) in the domain,
# and coloring does therefore not prevent race conditions during multithreaded assembly.
# For a domain buffer set up with `threading = true`, use
# `KeReAssembler(K, r; atomic = true)` to assemble domains with algebraic variables multithreaded.
# Otherwise, such domains are assembled sequentially.
