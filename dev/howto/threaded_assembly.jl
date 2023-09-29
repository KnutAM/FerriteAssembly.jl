using Ferrite, FerriteAssembly
import FerriteAssembly.ExampleElements: StationaryFourier

grid = generate_grid(Quadrilateral, (20, 20))
ip = Lagrange{RefQuadrilateral,1}()
dh = DofHandler(grid); add!(dh, :u, ip); close!(dh)
cellvalues = CellValues(QuadratureRule{RefQuadrilateral}(2), ip);
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
