# # Robin boundary conditions
# Robin boundary condition is a special type of boundary condition that may be interpreted as 
# a Neumann boundary condition, but where the flux depends on the field-variable. For scalar 
# problems, we can write the Robin BC contribution as 
# ```math
# \int_{\Gamma_\mathrm{R}} \delta u\ q_\mathrm{n}(u)\ \mathrm{d}\Gamma
# ```
# where $u$ is the primary variable, $\Gamma_\mathrm{R}$ is the Robin boundary, and 
# $q_\mathrm{n}(u)$ is the boundary flux that depends on the solution $u$. In this "how-to",
# we implement a simple linear Robin boundary condition, $q_\mathrm{n}(u) = k[u - u_\mathrm{b}]$.
# This can be interpreted physically for e.g. a thermal problem as a resistance, $1/k$, of heat transfer 
# from the domain to the outside, where the latter has a constant temperature $u_\mathrm{b}$. 

# ## Implementing Robin BC
# To implement this, we will overload the `face_routine!` function, which allows us to calculate 
# both a residual and stiffness contribution from faces. Hence, we first define a "material", 
# `RobinBC`, and implement the `face_routine!` that encodes the contribution to the weak form as 
# outlined above.  
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

# ## Standard `Ferrite.jl` setup
# To assemble the contributions from the Robin boundary, we need the regular setup
# for a Ferrite simulation, e.g. 
grid = generate_grid(Quadrilateral, (10,10))
ip = Lagrange{RefQuadrilateral,1}()

dh = DofHandler(grid)
add!(dh, :u, ip)
close!(dh)

K = create_sparsity_pattern(dh)
r = zeros(ndofs(dh))
a = zeros(ndofs(dh));

# ## Using Robin BC
# For the Robin boundary condition, we also need to define the integration via 
# FaceValues, in addition to an instance of `RobinBC`, for which we set the 
# "outside temperature" to -1.0. 
face_qr = FaceQuadratureRule{RefQuadrilateral}(2); 
fv = FaceValues(face_qr, ip);
rbc = RobinBC(1.0, -1.0);

# Finally, we set up a domain buffer with the actual materials and faceset,
# and create an assembler. Naturally, the boundary conditions are combined with
# other contributions, such as the assembly over cells. Typically, the same 
# assembler can be used. If not, it is important to set the keyword arguments 
# such that `K` and `r` are not zeroed between different calls to `work!`.
domainbuffer = setup_domainbuffer(DomainSpec(dh, rbc, fv; set=getfaceset(grid, "right")))
assembler = start_assemble(K, r)
work!(assembler, domainbuffer; a=a);

# ## Visualizing the BC
# To visualize the output, we can export the vtk the residual, showing how 
# contributions have been added to the boundary. 
vtk_grid("RobinBC", grid) do vtk
    vtk_point_data(vtk, dh, r)
end;

