using Ferrite
using FerriteAssembly
import FerriteAssembly.ExampleElements: StationaryFourier


function setup()
    grid = generate_grid(Quadrilateral, (2000,2000))
    ip = Ferrite.default_interpolation(getcelltype(grid))
    qr = QuadratureRule{2,RefCube}(2)
    dh = DofHandler(grid); add!(dh, :u, 1, ip); close!(dh)
    cv = CellScalarValues(qr, ip)

    buffer, new_states = setup_assembly(dh, StationaryFourier(1.0), cv; threading=Val(true))
    K = create_sparsity_pattern(dh)
    r = zeros(ndofs(dh))
    
    return K, r, new_states, buffer
end

const K_, r_, new_states_, buffer_ = setup()
const assembler_ = start_assemble(K, r)
@time doassemble!(assembler_, new_states_, buffer_)
@time doassemble!(assembler_, new_states_, buffer_)