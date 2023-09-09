# # Infinite plate with hole
#
# ![](plate_with_hole.png)
#-
# In this example we will solve a simple elasticity problem; an infinite plate with a hole.
# The main goal of the tutorial is to show how one can solve the problem using Isogeometric Analysis (IGA), 
# or in other words, solving a FE-problem with splines as the basis/shape functions.
# By using so called bezier extraction, we will see that most of the structure of the code will be the same as in standard FE-codes (however many differences are happening "under the hood").

#md # !!! note 
#md #     It is expected that the reader already be familiar with IGA and the concept of "bezier extraction".
#md #     It is also expected that the reader is familiar with the Ferrite package. In particular Ferrite.DofHandler and Ferrite.CellValues.

# Start by loading the necessary packages
using Ferrite, IGA, LinearAlgebra, FerriteAssembly
import FerriteAssembly.ExampleElements: ElasticPlaneStrain
using Serialization

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
    ## Performing the computation on a NURBS-patch is possible, 
    ## but it is much easier to use the bezier-extraction technique. 
    ## For this we transform the NURBS-patch into a `BezierGrid`. 
    ## The `BezierGrid` is identical to the standard `Ferrite.Grid`, 
    ## but includes the NURBS-weights and bezier extraction operators.
    grid = BezierGrid(nurbsmesh)

    ## For creating face sets, note that the nodes/controlpoints may 
    ## not be exactly on the geometry due to the non-interpolatory 
    ## nature of NURBS spline functions. However, in most cases 
    ## they will be close enough to use the Ferrite functions below.
    addnodeset!(grid,"right", (x) -> x[1] ≈ -0.0)
    addfaceset!(grid, "left", (x) -> x[1] ≈ -4.0)
    addfaceset!(grid, "bot", (x) -> x[2] ≈ 0.0)
    addfaceset!(grid, "right", (x) -> x[1] ≈ 0.0)

    return grid 
end

grid = create_mesh();

# ### Define interpolations, integration and FE-values
# Create the cellvalues storing the shape function values. 
# Note that the `CellVectorValues`/`FaceVectorValues` are wrapped in a `BezierValues`. 
# It is in the reinit-function of the `BezierValues` that the actual bezier transformation 
# of the shape values is performed. 
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
ch = ConstraintHandler(dh);
# The bottom face should only be able to move in x-direction
dbc1 = Dirichlet(:u, getfaceset(grid, "bot"), Returns(0.0), 2)
add!(ch, dbc1);

# The right boundary should only be able to move in y-direction
dbc2 = Dirichlet(:u, getfaceset(grid, "right"), Returns(0.0), 1)
add!(ch, dbc2)

close!(ch)
update!(ch, 0.0);

# Then Neumann conditions:
# Apply outwards traction on the left surface 
lh = LoadHandler(dh)
traction = Vec((-10.0, 0.0))
nbc = Neumann(:u, fv, getfaceset(grid, "left"), Returns(-traction))
add!(lh, nbc);

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
# Now we want to export the results to VTK. So we calculate the stresses in each gauss-point, and project them to 
# the nodes using the L2Projector from Ferrite. Node that we need to create new CellValues of type CellScalarValues, since the 
# L2Projector only works with scalar fields.
# Prepare postprocessing - i.e. calculating the stresses given the solution. 
struct CellStress{TT}
    s::Vector{Vector{TT}}
end
function CellStress(grid::Ferrite.AbstractGrid)
    Tb = SymmetricTensor{2,Ferrite.getdim(grid)}
    TT = Tb{Float64,Tensors.n_components(Tb)}
    return CellStress([TT[] for _ in 1:getncells(grid)])
end

function FerriteAssembly.integrate_cell!(cellstress::CellStress, state, ae, material, cv, cellbuffer)
    σ = cellstress.s[cellid(cellbuffer)]
    for q_point in 1:getnquadpoints(cv)
        ϵ = function_symmetric_gradient(cv, q_point, ae)
        push!(σ, material.C⊡ϵ)
    end
end

cellstresses = CellStress(grid)
integrator = Integrator(cellstresses)
work!(integrator, buffer; a=a)

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
