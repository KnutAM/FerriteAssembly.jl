using Ferrite, FerriteAssembly
import FerriteAssembly.ExampleElements: ElasticPlaneStrain
ip = Lagrange{RefQuadrilateral,1}()^2
grid = generate_grid(Quadrilateral, (3,3))

dh = DofHandler(grid)
add!(dh, :u, ip)
close!(dh)

ch = ConstraintHandler(dh)
add!(ch, Dirichlet(:u, getfacetset(grid, "left"), Returns(zero(Vec{2}))))
add!(ch, Dirichlet(:u, getfacetset(grid, "right"), Returns(1.0), [1,]))
close!(ch)
update!(ch, 0.0)

qr = QuadratureRule{RefQuadrilateral}(2)
cv = CellValues(qr, ip, ip)
m = ElasticPlaneStrain(;E=80e3, ν=0.3)

buffer = setup_domainbuffer(DomainSpec(dh, m, cv))
K = create_sparsity_pattern(dh)
r = zeros(ndofs(dh))
a = zeros(ndofs(dh));

apply!(a, ch)
assembler = KeReAssembler(K, r; ch, apply_zero=true)
work!(assembler, buffer; a=a);

a .-= K\r;

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
