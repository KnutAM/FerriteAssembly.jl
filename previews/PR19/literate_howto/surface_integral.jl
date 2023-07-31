# # Surface integration
using Ferrite, FerriteAssembly
# This how-to shows integration of the normal flux 
# on a surface. As usual, we need the basic Ferrite 
# building blocks, in this case a grid, dofhandler, 
# and facevalues 
grid = generate_grid(Hexahedron, (5,5,5))
dh = DofHandler(grid); add!(dh, :u, 1); close!(dh)
qr = QuadratureRule{2,RefCube}(2); ip = Lagrange{3,RefCube,1}()
fv = FaceScalarValues(qr, ip);

# We also need a solution vector to integrate, 
# unless we only calculate geometric properties.
a = zeros(ndofs(dh))
apply_analytical!(a, dh, :u, norm) # f(x)=norm(x)

# And then we decide which faces to integrate over,
# and define the integrator, here using the simple
# interface 
domainbuffer = setup_domainbuffer(DomainSpec(dh, nothing, fv; set=getfaceset(grid, "right")))
integrator = SimpleIntegrator((u,∇u,n)->∇u⋅n, 0.0)
work!(integrator, domainbuffer; a=a)

println("Flux: ∫qₙ dA = ", integrator.val)