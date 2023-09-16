using Ferrite, FerriteAssembly
import FerriteAssembly.ExampleElements: StationaryFourier

grid = generate_grid(Quadrilateral, (20, 20))
dh = DofHandler(grid); add!(dh, :u, 1); close!(dh)
cellvalues = CellScalarValues(QuadratureRule{2, RefCube}(2), Lagrange{2, RefCube, 1}());
K = create_sparsity_pattern(dh)
r = zeros(ndofs(dh));

material = StationaryFourier(1.0);

domain = DomainSpec(dh, material, cellvalues)
buffer = setup_domainbuffer(domain; threading=true);

assembler = start_assemble(K, r);

work!(assembler, buffer);

lh = LoadHandler(dh)
add!(lh, BodyLoad(:u, 2, Returns(-1.0))) # rᵢ -= ∫ δuᵢ*1.0*dV
apply!(r, lh, 0.0) # t = 0

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
