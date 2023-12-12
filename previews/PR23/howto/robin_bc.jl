using Ferrite, FerriteAssembly

struct RobinBC{T}
    k::T
    ub::T
end

function FerriteAssembly.face_routine!(
        Ke, re, ae, rbc::RobinBC, fv::FaceValues, facebuffer
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

K = create_sparsity_pattern(dh)
r = zeros(ndofs(dh))
a = zeros(ndofs(dh));

face_qr = FaceQuadratureRule{RefQuadrilateral}(2);
fv = FaceValues(face_qr, ip);
rbc = RobinBC(1.0, -1.0);

domainbuffer = setup_domainbuffer(DomainSpec(dh, rbc, fv; set=getfaceset(grid, "right")))
assembler = start_assemble(K, r)
work!(assembler, domainbuffer; a=a);

vtk_grid("RobinBC", grid) do vtk
    vtk_point_data(vtk, dh, r)
end;

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
