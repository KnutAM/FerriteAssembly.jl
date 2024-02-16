using Ferrite
using FerriteAssembly
import FerriteAssembly.ExampleElements: LinearElastic

threading = Threads.nthreads() > 0
n = 30
grid = generate_grid(Hexahedron, (10*n,n,n))

ip = Lagrange{RefHexahedron,1}()^3
qr = QuadratureRule{RefHexahedron}(2)
dh = DofHandler(grid)
add!(dh, :u, ip)
close!(dh)
cv = CellValues(qr, ip)
material = LinearElastic(E=200e9, ν=0.3)
buffer = setup_domainbuffer(DomainSpec(dh, material, cv); threading)

K = create_sparsity_pattern(dh)
r = zeros(ndofs(dh))
assembler = start_assemble(K, r)

work!(assembler, buffer)

function dotiming(K, r, buffer)
    assembler = start_assemble(K, r)
    @elapsed work!(assembler, buffer)
end

open("bench_thread_ferrite_assembly.txt"; append=true) do fid
    write(fid, string(Threads.nthreads()))
    write(fid, " ")
    for _ in 1:3
        t = dotiming(K, r, buffer)
        write(fid, string(t))
        write(fid, " ")
    end
    write(fid, "\n")
end
