# # Local constraint application
# **Executive summary:** Replace `start_assemble(K,r)` with `KeReAssembler(K,r;ch,apply_zero=true)`.
# 
# In some cases, it might be beneficial to apply constraints locally 
# by using `Ferrite`'s `apply_assemble!`. This is supported by using 
# `FerriteAssembly`'s `ReAssembler` and `KeReAssembler`.
# To demonstrate, let's start by setting up a quick simulation setup
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

# We can now create the special `KeReAssembler`, 
# which also accepts the constraint handler as input, 
# in addition to the `apply_zero` keyword that is 
# forwarded to `Ferrite.apply_assemble!`.
apply!(a, ch)
assembler = KeReAssembler(K, r; ch, apply_zero=true)
work!(assembler, buffer; a=a);

# And finally we can solve our problem update 
a .-= K\r;

using Test                              #src
K2 = similar(K)                         #src
r2 = similar(r)                         #src
a2 = zeros(ndofs(dh))                   #src
apply!(a2, ch)                          #src
std_assembler = start_assemble(K2, r2)  #src
work!(std_assembler, buffer; a=a2)      #src
apply!(K2, r2, ch)                      #src
a2 .-= K\r                              #src
@test a2 ≈ a                            #src
