using Ferrite, FerriteAssembly

grid = generate_grid(Hexahedron, (5,5,5))
ip = Lagrange{RefHexahedron,1}()
dh = DofHandler(grid); add!(dh, :u, ip); close!(dh)
qr = FacetQuadratureRule{RefHexahedron}(2);
fv = FacetValues(qr, ip);

a = zeros(ndofs(dh))
apply_analytical!(a, dh, :u, norm); # f(x)=norm(x)

domainbuffer = setup_domainbuffer(DomainSpec(dh, nothing, fv; set=getfacetset(grid, "right")));

s_integrator = SimpleIntegrator((u,∇u,n)->∇u⋅n, 0.0);

work!(s_integrator, domainbuffer; a=a);

println("Flux: ∫qₙ dA = ", s_integrator.val)

mutable struct NormalFlux{T}
    qn::T
end
function FerriteAssembly.integrate_facet!(nf::NormalFlux, ae, material, fv, facetbuffer)
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
