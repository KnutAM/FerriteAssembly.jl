using Ferrite
using FerriteAssembly

# Implementation of LinearStiffness material which only calculates stiffness and external foreach
# due to a volumetric force vector `b` stored in `user_data`
struct LinearStiffness{TT<:SymmetricTensor{4}}
    C::TT
end
function LinearStiffness(;E, ν)
    λ = E*ν / ((1+ν) * (1 - 2ν))
    μ = E / (2(1+ν))
    δ(i,j) = i == j ? 1.0 : 0.0
    g(i,j,k,l) = λ*δ(i,j)*δ(k,l) + μ*(δ(i,k)*δ(j,l) + δ(i,l)*δ(j,k))
    C = SymmetricTensor{4, 3}(g);
    return LinearStiffness(C)
end
FerriteAssembly.allocate_cell_cache(::LinearStiffness, cv) = [shape_symmetric_gradient(cv, 1, i) for i in 1:getnbasefunctions(cv)]

function FerriteAssembly.element_routine!(Ke, re, state, ae, m::LinearStiffness, cv, buffer)
    ɛ = FerriteAssembly.get_user_cache(buffer)
    b = FerriteAssembly.get_user_data(buffer)
    n_basefuncs = getnbasefunctions(cv)
    for q_point in 1:getnquadpoints(cv)
        for i in 1:n_basefuncs
            ɛ[i] = symmetric(shape_gradient(cv, q_point, i))
        end
        dΩ = getdetJdV(cv, q_point)
        for i in 1:n_basefuncs
            δu = shape_value(cv, q_point, i)
            re[i] -= (δu ⋅ b) * dΩ
            ɛC = ɛ[i] ⊡ m.C
            for j in 1:n_basefuncs
                Ke[i, j] += (ɛC ⊡ ɛ[j]) * dΩ
            end
        end
    end
end
# End of LinearStiffness implementation

threading = Threads.nthreads() > 0
n = 50
grid = generate_grid(Hexahedron, (10*n,n,n))

ip = Lagrange{RefHexahedron,1}()^3
qr = QuadratureRule{RefHexahedron}(2)
dh = DofHandler(grid)
add!(dh, :u, ip)
close!(dh)
cv = CellValues(qr, ip)
material = LinearStiffness(; E=200e9, ν=0.3)
b = Vec{3}((0.0, 0.0, 0.0))
buffer = setup_domainbuffer(DomainSpec(dh, material, cv; user_data = b); threading)

K = create_sparsity_pattern(dh)
r = zeros(ndofs(dh))
assembler = start_assemble(K, r)

work!(assembler, buffer)

function dotiming(K, r, buffer)
    assembler = start_assemble(K, r)
    @elapsed work!(assembler, buffer)
end

open("threaded_elasticity2.txt"; append=true) do fid
    write(fid, string(Threads.nthreads()))
    write(fid, " ")
    for _ in 1:3
        t = dotiming(K, r, buffer)
        write(fid, string(t))
        write(fid, " ")
    end
    write(fid, "\n")
end
