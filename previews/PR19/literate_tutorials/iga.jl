# # Isogeometric analysis with `IGA.jl`
# This tutorial shows how to integrate FerriteAssembly with the 
# isogeometric analysis toolbox IGA.jl, directly based on `IGA.jl`'s 
# [Infinite plate with hole](https://lijas.github.io/IGA.jl/dev/examples/plate_with_hole/)
# example. Hence, please see there for further documentation details and important remarks 
# regarding IGA. 
# 
# The example considers solving a plate with a hole. A quater of a plate is considered via symmetry 
# boundary conditions, and a tensile load is applied on one edge.

# Start by loading the necessary packages
using Ferrite, IGA, LinearAlgebra, FerriteAssembly
import FerriteAssembly.ExampleElements: ElasticPlaneStrain

# ## Extensions and fixes to IGA.jl 
Ferrite.getnormal(::BezierFaceValues{sdim,T}, ::Int) where {sdim,T} = T(NaN)*zero(Vec{sdim,T})
function Ferrite.spatial_coordinate(bv::IGA.BezierValues, q_point::Int, bc::BezierCoords)
    return spatial_coordinate(bv, q_point, (bc.xb, bc.wb))
end
function Ferrite.function_symmetric_gradient(bv::IGA.BezierValues, q_point::Int, u::AbstractVector, args...)
    return function_symmetric_gradient(bv.cv_store, q_point, u, args...)
end
function Ferrite.shape_symmetric_gradient(bv::IGA.BezierValues, q_point::Int, i::Int)
    return shape_symmetric_gradient(bv.cv_store, q_point, i)
end

# This looks like a bug, and need to just overload a more specific version.
function Ferrite.getcoordinates(grid::BezierGrid{dim,C,T}, ic::Int) where {dim,C,T<:Float64}
    n = Ferrite.nnodes_per_cell(grid, ic)
    w = zeros(T, n)
    wb = zeros(T, n)
    xb = zeros(Vec{dim,T}, n)
    x = zeros(Vec{dim,T}, n)
    
    ## `copy` to avoid reference the extraction operator
    bc = BezierCoords(xb, wb, x, w, copy(grid.beo[ic]))

    return getcoordinates!(bc,grid,ic)
end
function Ferrite.getcoordinates!(bc::BezierCoords{dim,Float64}, grid::BezierGrid, ic::Int) where {dim}
    get_bezier_coordinates!(bc.xb, bc.wb, bc.x, bc.w, grid, ic)
    copyto!(bc.beo, get_extraction_operator(grid, ic)) 
    return bc
end

# ## Problem setup
# ### Mesh
# We begin by generating the mesh by using `IGA.jl`'s built-in mesh generators 
# In this example, we will generate the patch called "plate with hole". 
# Note, currently this function can only generate the patch with second order basefunctions. 
function create_mesh(;nels=(20,10))
    nurbsmesh = generate_nurbs_patch(:plate_with_hole, nels) 
    grid = BezierGrid(nurbsmesh)

    addfaceset!(grid, "left", (x) -> x[1] ≈ -4.0)
    addfaceset!(grid, "bot", (x) -> x[2] ≈ 0.0)
    addfaceset!(grid, "right", (x) -> x[1] ≈ 0.0)

    return grid 
end
grid = create_mesh();

# ### Define interpolations, integration and FE-values
# Create the cellvalues storing the shape function values. 
orders=(2,2)
ip = BernsteinBasis{2,orders}()
qr_cell = QuadratureRule{2,RefCube}(4)
qr_face = QuadratureRule{1,RefCube}(3)

cv = BezierCellValues( CellVectorValues(qr_cell, ip) );
fv = BezierFaceValues( FaceVectorValues(qr_face, ip) );

# ### Dof distribution
# Distribute dofs as normal
dh = MixedDofHandler(grid)
push!(dh, :u, 2, ip)
close!(dh);

# And allocate system matrices and vectors
a = zeros(ndofs(dh))
r = zeros(ndofs(dh))
K = create_sparsity_pattern(dh);

# ### Boundary conditions 
# Starting with Dirichlet conditions: 
# 1) Bottom face should only be able to move in x-direction
# 2) Right boundary should only be able to move in y-direction
ch = ConstraintHandler(dh);
add!(ch, Dirichlet(:u, getfaceset(grid, "bot"), Returns(0.0), 2))
add!(ch, Dirichlet(:u, getfaceset(grid, "right"), Returns(0.0), 1))
close!(ch)
update!(ch, 0.0);

# Then Neumann conditions:
# Apply outwards traction on the left surface,
# taking the negative value since r = fint - fext.
traction = Vec((-10.0, 0.0))
lh = LoadHandler(dh)
add!(lh, Neumann(:u, fv, getfaceset(grid, "left"), Returns(-traction)));

# ### FerriteAssembly setup 
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

# using Serialization #src
# uref = deserialize("uref.bin") #src
# sref = deserialize("sref.bin") #src
# @show uref ≈ a #src
# @show sref ≈ σ_nodes #src
