# # Surface integration
using Ferrite, FerriteAssembly
# This how-to shows integration of the normal flux 
# on a surface. As usual, we need the basic Ferrite 
# building blocks, in this case a grid, dofhandler, 
# and facetvalues 
grid = generate_grid(Hexahedron, (5,5,5))
ip = Lagrange{RefHexahedron,1}()
dh = DofHandler(grid); add!(dh, :u, ip); close!(dh)
qr = FacetQuadratureRule{RefHexahedron}(2); 
fv = FacetValues(qr, ip);

# We also need a solution vector to integrate, 
# unless we only calculate geometric properties.
a = zeros(ndofs(dh))
apply_analytical!(a, dh, :u, norm); # f(x)=norm(x)

# And then we decide which facets to integrate over
domainbuffer = setup_domainbuffer(DomainSpec(dh, nothing, fv; set=getfacetset(grid, "right")));

# ## Using `SimpleIntegrator`
# Using the simple interface, `SimpleIntegrator`, we simply
# give it the function, `(u,∇u,n)->∇u⋅n`, and the initial value
s_integrator = SimpleIntegrator((u,∇u,n)->∇u⋅n, 0.0);

# And then let it do its work
work!(s_integrator, domainbuffer; a=a);

println("Flux: ∫qₙ dA = ", s_integrator.val)

# ## Using `Integrator`
# To demonstrate how the full-fledged interface can be used, we 
# perform the same task by using the `Integrator`.
# To do that, we have to define an integrand, let's call it `NormalFlux`,
# and overload the facet integration function
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

# To do the actual integration, we define an instance of the integrand,
# create an `Integrator`, and do the work.
nf = NormalFlux(0.0)
integrator = Integrator(nf)
work!(integrator, domainbuffer; a=a)

println("Flux: ∫qₙ dA = ", nf.qn)