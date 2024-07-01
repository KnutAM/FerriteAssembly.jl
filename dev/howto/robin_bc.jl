using Ferrite, FerriteAssembly

struct RobinBC{T}
    k::T
    ub::T
end

function FerriteAssembly.facet_routine!(
        Ke, re, ae, rbc::RobinBC, fv::FacetValues, facetbuffer
        )
    for q_point in 1:getnquadpoints(fv)
        dΓ = getdetJdV(fv, q_point)
        u = function_value(fv, q_point, ae)
        qn = rbc.k*(u - rbc.ub) # Flux out of domain
        for i in 1:getnbasefunctions(fv)
            δN = shape_value(fv, q_point, i)
            re[i] += δN*qn*dΓ
            for j in 1:getnbasefunctions(fv)
                N = shape_value(fv, q_point, j)
                Ke[i,j] += δN*rbc.k*N*dΓ
            end
        end
    end
end;

grid = generate_grid(Quadrilateral, (10,10))
ip = Lagrange{RefQuadrilateral,1}()

dh = DofHandler(grid)
add!(dh, :u, ip)
close!(dh)

K = allocate_matrix(dh)
r = zeros(ndofs(dh))
a = zeros(ndofs(dh));

facet_qr = FacetQuadratureRule{RefQuadrilateral}(2);
fv = FacetValues(facet_qr, ip);
rbc = RobinBC(1.0, -1.0);

domainbuffer = setup_domainbuffer(DomainSpec(dh, rbc, fv; set=getfacetset(grid, "right")))
assembler = start_assemble(K, r)
work!(assembler, domainbuffer; a=a);

VTKFile("RobinBC", grid) do vtk
    write_solution(vtk, dh, r)
end;

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
