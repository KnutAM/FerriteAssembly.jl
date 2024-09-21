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

function calculate_stress(m, u, ∇u, qp_state)
    ϵ = symmetric(∇u)
    return m.C ⊡ ϵ
end

qe = QuadPointEvaluator{SymmetricTensor{2,2,Float64,3}}(buffer, calculate_stress);
work!(qe, buffer; a=a)

IGA.VTKIGAFile("plate_with_hole.vtu", grid) do vtk
    write_solution(vtk, dh, a)
end

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
