using Ferrite, IGA, LinearAlgebra, FerriteAssembly
import FerriteAssembly.ExampleElements: ElasticPlaneStrain

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

ip = IGAInterpolation{RefQuadrilateral, order}()
qr_cell = QuadratureRule{RefQuadrilateral}(4)
qr_face = FacetQuadratureRule{RefQuadrilateral}(3)

cv = BezierCellValues(qr_cell, ip^2);
fv = BezierFacetValues(qr_face, ip^2);

dh = DofHandler(grid)
add!(dh, :u, ip^2)
close!(dh);

a = zeros(ndofs(dh))
r = zeros(ndofs(dh))
K = allocate_matrix(dh);

ch = ConstraintHandler(dh);
add!(ch, Dirichlet(:u, getfacetset(grid, "bot"), Returns(0.0), 2))
add!(ch, Dirichlet(:u, getfacetset(grid, "right"), Returns(0.0), 1))
close!(ch)
update!(ch, 0.0);

traction = Vec((-10.0, 0.0))
lh = LoadHandler(dh)
add!(lh, Neumann(:u, fv, getfacetset(grid, "left"), Returns(-traction)));

material = ElasticPlaneStrain(;E=100.0, ν=0.3)
domain = DomainSpec(dh, material, cv)
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
    TT = typeof(zero(SymmetricTensor{2, Ferrite.getspatialdim(grid)}))
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

IGA.VTKIGAFile("plate_with_hole.vtu", grid) do vtk
    write_solution(vtk, dh, a)
end

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
