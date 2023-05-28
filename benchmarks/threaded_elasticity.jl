using Ferrite
using FerriteAssembly
import FerriteAssembly.ExampleElements: LinearElastic

setup(threaded=Val(Threads.nthreads()>1)) = setup_threaded(threaded)
function setup_threaded(threaded::Val)
    grid = generate_grid(Hexahedron, (20,20,20))
    ip = Ferrite.default_interpolation(getcelltype(grid))
    qr = QuadratureRule{3,RefCube}(5)
    dh = DofHandler(grid); add!(dh, :u, 3, ip); close!(dh)
    cv = CellVectorValues(qr, ip)

    buffer, new_states = setup_assembly(dh, LinearElastic(E=210e9, ν=0.3), cv; threading=threaded)
    K = create_sparsity_pattern(dh)
    r = zeros(ndofs(dh))
    
    return K, r, new_states, buffer
end

K, r, new_states, buffer = setup()
assembler = start_assemble(K, r)
doassemble!(assembler, new_states, buffer)

#assembler = start_assemble(K, r)
#@time doassemble!(assembler, new_states, buffer)

function dotiming(K, r, new_states, buffer)
    assembler = start_assemble(K, r)
    @time doassemble!(assembler, new_states, buffer)
end

for i in 1:5
    print(i, ": "); dotiming(K, r, new_states, buffer)
end