using Ferrite, FerriteAssembly

grid = generate_grid(Hexahedron, (5,5,5))
dh = DofHandler(grid); add!(dh, :u, 1); close!(dh)
qr = QuadratureRule{2,RefCube}(2); ip = Lagrange{3,RefCube,1}()
fv = FaceScalarValues(qr, ip);

a = zeros(ndofs(dh))
apply_analytical!(a, dh, :u, norm) # f(x)=norm(x)

domainbuffer = setup_domainbuffer(DomainSpec(dh, nothing, fv; set=getfaceset(grid, "right")))
integrator = SimpleIntegrator((u,∇u,n)->∇u⋅n, 0.0)
work!(integrator, domainbuffer; a=a)

println("Flux: ∫qₙ dA = ", integrator.val)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

