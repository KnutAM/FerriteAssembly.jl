using Ferrite, FerriteAssembly

grid = generate_grid(Hexahedron, (5,5,5))
dh = DofHandler(grid); add!(dh, :u, 1); close!(dh)
qr = QuadratureRule{2,RefCube}(2); ip = Lagrange{3,RefCube,1}()
fv = FaceScalarValues(qr, ip);

a = zeros(ndofs(dh))
apply_analytical!(a, dh, :u, norm); # f(x)=norm(x)

domainbuffer = setup_domainbuffer(DomainSpec(dh, nothing, fv; set=getfaceset(grid, "right")));

s_integrator = SimpleIntegrator((u,∇u,n)->∇u⋅n, 0.0);

work!(s_integrator, domainbuffer; a=a);

println("Flux: ∫qₙ dA = ", s_integrator.val)

mutable struct NormalFlux{T}
    qn::T
end
function FerriteAssembly.integrate_face!(nf::NormalFlux, ae, material, fv, facebuffer)
    for q_point in 1:getnquadpoints(fv)
        dA = getdetJdV(fv, q_point)
        ∇u = function_gradient(fv, q_point, ae)
        n = getnormal(fv, q_point)
        nf.qn += (∇u⋅n)*dA
    end
end;

nf = NormalFlux(0.0)
integrator = Integrator(nf)
work!(integrator, domainbuffer; a=a)

println("Flux: ∫qₙ dA = ", nf.qn)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
