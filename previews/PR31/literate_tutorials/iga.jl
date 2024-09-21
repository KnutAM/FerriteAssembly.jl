# # Using `IGA.jl`
# This tutorial shows how to integrate FerriteAssembly with the 
# isogeometric analysis toolbox IGA.jl, directly based on `IGA.jl`'s 
# [Infinite plate with hole](https://lijas.github.io/IGA.jl/dev/examples/plate_with_hole/)
# example. Hence, please see there for further documentation details and important remarks 
# regarding IGA. 
# 
# The example considers solving a plate with a hole. A quarter of a plate is considered via symmetry 
# boundary conditions, and a tensile load is applied on one edge.
# The full script without intermediate comments is available at the 
# [bottom of this page](@ref iga_plain_program).

# Start by loading the necessary packages
using Ferrite, IGA, LinearAlgebra, FerriteAssembly
import FerriteAssembly.ExampleElements: ElasticPlaneStrain

# ## Setup
# To clarify the differences, we split the setup into `IGA.jl`, `Ferrite.jl`,
# and `FerriteAssembly.jl` specific parts.
# ### `IGA.jl` setup 
# We begin by generating the mesh by using `IGA.jl`'s built-in mesh generators, 
# specifically a "plate with hole". 
function create_mesh(; nels = (20,10), order)
    nurbsmesh = generate_nurbs_patch(:plate_with_hole, nels, order)
    grid = BezierGrid(nurbsmesh)

    addfacetset!(grid, "left", (x) -> x[1] ≈ -4.0)
    addfacetset!(grid, "bot", (x) -> x[2] ≈ 0.0)
    addfacetset!(grid, "right", (x) -> x[1] ≈ 0.0)

    return grid 
end
order = 2
grid = create_mesh(; order);

# Create the `IGA.jl` cell and facet values for storing the 
# the `IGA.jl` shape function values and gradients. 
ip = IGAInterpolation{RefQuadrilateral, order}()
qr_cell = QuadratureRule{RefQuadrilateral}(4)
qr_face = FacetQuadratureRule{RefQuadrilateral}(3)

cv = BezierCellValues(qr_cell, ip^2);
fv = BezierFacetValues(qr_face, ip^2);

# ### `Ferrite.jl` setup
# Distribute dofs as normal
dh = DofHandler(grid)
add!(dh, :u, ip^2)
close!(dh);

# And allocate system matrices and vectors
a = zeros(ndofs(dh))
r = zeros(ndofs(dh))
K = allocate_matrix(dh);

# Before adding boundary conditions, starting with Dirichlet: 
# 1) Bottom face should only be able to move in x-direction
# 2) Right boundary should only be able to move in y-direction
ch = ConstraintHandler(dh);
add!(ch, Dirichlet(:u, getfacetset(grid, "bot"), Returns(0.0), 2))
add!(ch, Dirichlet(:u, getfacetset(grid, "right"), Returns(0.0), 1))
close!(ch)
update!(ch, 0.0);

# ### `FerriteAssembly.jl` setup 
# Neumann boundary conditions are added using `FerriteAssembly`'s 
# `LoadHandler`. We apply outwards traction on the left surface,
# and take the negative value since r = fint - fext.
traction = Vec((-10.0, 0.0))
lh = LoadHandler(dh)
add!(lh, Neumann(:u, fv, getfacetset(grid, "left"), Returns(-traction)));

material = ElasticPlaneStrain(;E=100.0, ν=0.3)
domain = DomainSpec(dh, material, cv)
buffer = setup_domainbuffer(domain);

# ## Assembly and solution
# We first assemble the equation system,
assembler = start_assemble(K, r)
apply!(a, ch)
work!(assembler, buffer; a=a)
apply!(r, lh, 0.0);

# before solving it,
apply_zero!(K, r, ch)
a .-= K\r
apply!(a, ch);

# ## Postprocessing
# We use the `QuadPointEvaluator` to calculate stresses in each integration point,
# The `QuadPointEvaluator` requires a function with the following input arguments
function calculate_stress(m, u, ∇u, qp_state)
    ϵ = symmetric(∇u)
    return m.C ⊡ ϵ
end;

# We can then use this to calculate the stress tensor in each integration point
qe = QuadPointEvaluator{SymmetricTensor{2,2,Float64,3}}(buffer, calculate_stress);
work!(qe, buffer; a=a);

# Now we want to export the results to VTK. So we project the stresses at 
# the quadrature points to the nodes using the L2Projector from Ferrite. 
# Currently, however, IGA doesn't support L2 projection. 
# ```julia
# # projector = L2Projector(ip, grid)
# # σ_nodes = IGA.igaproject(projector, qe.data, qr_cell; project_to_nodes=true);
# ```

# Output results to VTK
IGA.VTKIGAFile("plate_with_hole.vtu", grid) do vtk
    write_solution(vtk, dh, a)
end

using Test                                      #src
# @test sum(norm, σ_nodes) ≈ 3087.2447327126742 #src
@test norm(norm.(qe.data)) ≈ 679.3207411544098  #src

#md # ## [Plain program](@id iga_plain_program)
#md #
#md # Here follows a version of the program without any comments.
#md # The file is also available here: [`iga.jl`](iga.jl).
#md #
#md # ```julia
#md # @__CODE__
#md # ```
