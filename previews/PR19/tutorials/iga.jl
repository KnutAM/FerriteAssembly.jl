using Ferrite, IGA, LinearAlgebra, FerriteAssembly
import FerriteAssembly.ExampleElements: ElasticPlaneStrain

function create_mesh(;nels=(20,10))
    nurbsmesh = generate_nurbs_patch(:plate_with_hole, nels)
    grid = BezierGrid(nurbsmesh)

    addfaceset!(grid, "left", (x) -> x[1] ≈ -4.0)
    addfaceset!(grid, "bot", (x) -> x[2] ≈ 0.0)
    addfaceset!(grid, "right", (x) -> x[1] ≈ 0.0)

    return grid
end
grid = create_mesh();

orders=(2,2)
ip = BernsteinBasis{2,orders}()
qr_cell = QuadratureRule{2,RefCube}(4)
qr_face = QuadratureRule{1,RefCube}(3)

cv = BezierCellValues( CellVectorValues(qr_cell, ip) );
fv = BezierFaceValues( FaceVectorValues(qr_face, ip) );

dh = MixedDofHandler(grid)
push!(dh, :u, 2, ip)
close!(dh);

a = zeros(ndofs(dh))
r = zeros(ndofs(dh))
K = create_sparsity_pattern(dh);

ch = ConstraintHandler(dh);
add!(ch, Dirichlet(:u, getfaceset(grid, "bot"), Returns(0.0), 2))
add!(ch, Dirichlet(:u, getfaceset(grid, "right"), Returns(0.0), 1))
close!(ch)
update!(ch, 0.0);

traction = Vec((-10.0, 0.0))
lh = LoadHandler(dh)
add!(lh, Neumann(:u, fv, getfaceset(grid, "left"), Returns(-traction)));

material = ElasticPlaneStrain(;E=100.0, ν=0.3)
domain = DomainSpec(FerriteAssembly.SubDofHandler(dh, dh.fieldhandlers[1]), material, cv)
buffer = setup_domainbuffer(domain);

assembler = start_assemble(K, r)
apply!(a, ch)
work!(assembler, buffer; a=a)
apply!(r, lh, 0.0);

apply_zero!(K, r, ch)
a .-= K\r
apply!(a, ch);

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

projector = L2Projector(ip, grid)
σ_nodes = IGA.igaproject(projector, cellstresses.s, qr_cell; project_to_nodes=true);

vtkgrid = vtk_grid("plate_with_hole.vtu", grid)
vtk_point_data(vtkgrid, dh, a)
vtk_point_data(vtkgrid, σ_nodes, "sigma", grid)
vtk_save(vtkgrid);

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
