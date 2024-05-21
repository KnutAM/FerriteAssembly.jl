# # Isogeometric analysis with `IGA.jl`
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
# We begin by generating the mesh by using `IGA.jl`'s built-in mesh generators 
# In this example, we will generate the patch called "plate with hole". 
# Note, currently this function can only generate the patch with second order basefunctions. 
function create_mesh(;nels=(20,10))
    nurbsmesh = generate_nurbs_patch(:plate_with_hole, nels) 
    grid = BezierGrid(nurbsmesh)

    addfacetset!(grid, "left", (x) -> x[1] ≈ -4.0)
    addfacetset!(grid, "bot", (x) -> x[2] ≈ 0.0)
    addfacetset!(grid, "right", (x) -> x[1] ≈ 0.0)

    return grid 
end
grid = create_mesh();

# Create the cellvalues storing the shape function values. 
orders=(2,2)
ip = BernsteinBasis{2,orders}()
qr_cell = QuadratureRule{2,RefCube}(4)
qr_facet = QuadratureRule{1,RefCube}(3)

cv = BezierCellValues( CellValues(qr_cell, ip^2) );
fv = BezierFacetValues( FacetValues(qr_facet, ip^2) );

# Distribute dofs as normal
dh = DofHandler(grid)
add!(dh, :u, ip^2)
close!(dh);

# And allocate system matrices and vectors
a = zeros(ndofs(dh))
r = zeros(ndofs(dh))
K = create_sparsity_pattern(dh);

# Starting with Dirichlet conditions: 
# 1) Bottom facet should only be able to move in x-direction
# 2) Right boundary should only be able to move in y-direction
ch = ConstraintHandler(dh);
add!(ch, Dirichlet(:u, getfacetset(grid, "bot"), Returns(0.0), 2))
add!(ch, Dirichlet(:u, getfacetset(grid, "right"), Returns(0.0), 1))
close!(ch)
update!(ch, 0.0);

# Then Neumann conditions:
# Apply outwards traction on the left surface,
# taking the negative value since r = fint - fext.
traction = Vec((-10.0, 0.0))
lh = LoadHandler(dh)
add!(lh, Neumann(:u, fv, getfacetset(grid, "left"), Returns(-traction)));

# FerriteAssembly setup 
material = ElasticPlaneStrain(;E=100.0, ν=0.3)
domain = DomainSpec(FerriteAssembly.SubDofHandler(dh, dh.fieldhandlers[1]), material, cv)
buffer = setup_domainbuffer(domain);

# ## Assemble 
assembler = start_assemble(K, r)
apply!(a, ch)
work!(assembler, buffer; a=a)
apply!(r, lh, 0.0);

# ## Solve
apply_zero!(K, r, ch)
a .-= K\r
apply!(a, ch);

# ## Postprocessing
# First, the stresses in each integration point are calculated by using the Integrator.  
struct CellStress{TT}
    s::Vector{Vector{TT}}
end
function CellStress(grid::Ferrite.AbstractGrid)
    Tb = SymmetricTensor{2,Ferrite.getdim(grid)}
    TT = Tb{Float64,Tensors.n_components(Tb)}
    return CellStress([TT[] for _ in 1:getncells(grid)])
end

function FerriteAssembly.integrate_cell!(stress::CellStress, state, ae, material, cv, cb)
    σ = stress.s[cellid(cb)]
    for q_point in 1:getnquadpoints(cv)
        ϵ = function_symmetric_gradient(cv, q_point, ae)
        push!(σ, material.C⊡ϵ)
    end
end
cellstresses = CellStress(grid)
integrator = Integrator(cellstresses)
work!(integrator, buffer; a=a);

# Now we want to export the results to VTK. So we project the stresses at 
# the quadrature points to the nodes using the L2Projector from Ferrite. 
projector = L2Projector(ip, grid)
σ_nodes = IGA.igaproject(projector, cellstresses.s, qr_cell; project_to_nodes=true);

# Output results to VTK
vtkgrid = vtk_grid("plate_with_hole.vtu", grid)
vtk_point_data(vtkgrid, dh, a)
vtk_point_data(vtkgrid, σ_nodes, "sigma", grid)
vtk_save(vtkgrid);

using Test                                    #src
@test sum(norm, σ_nodes) ≈ 3087.2447327126742 #src

#md # ## [Plain program](@id iga_plain_program)
#md #
#md # Here follows a version of the program without any comments.
#md # The file is also available here: [`iga.jl`](iga.jl).
#md #
#md # ```julia
#md # @__CODE__
#md # ```
