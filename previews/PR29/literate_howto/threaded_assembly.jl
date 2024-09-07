# # Threaded assembly 
# **Executive summary:** Call `setup_domainbuffer(domain; threading=true)`.
# 
# FerriteAssembly aims to provide most functionality to be optionally done threaded, 
# to make sure that the assembly (or more generally, calls to `work!`) is not the bottleneck
# in the simulation time. In this example, we show how to do the assembly of the stationary 
# heat flow threaded, by using the example element, `StationaryFourier`: 
using Ferrite, FerriteAssembly
import FerriteAssembly.ExampleElements: StationaryFourier

# The setup is the same as always, i.e.
grid = generate_grid(Quadrilateral, (20, 20))
ip = Lagrange{RefQuadrilateral,1}()
dh = DofHandler(grid); add!(dh, :u, ip); close!(dh)
cellvalues = CellValues(QuadratureRule{RefQuadrilateral}(2), ip);
K = allocate_matrix(dh)
r = zeros(ndofs(dh));

# We then create an instance of the material, for which an element routine is already provided,
# but we could just as well have used the one form the heat equation tutorial. 
material = StationaryFourier(1.0);

# To set up a threaded assembly, we just need to set `threading=true` when setting up the domainbuffer.
domain = DomainSpec(dh, material, cellvalues)
buffer = setup_domainbuffer(domain; threading=true);
# This creates a default coloring of the grid, but custom coloring can also be given.

# The `worker` in this case is the standard Ferrite assembler:
assembler = start_assemble(K, r);

# Given this worker, we can do the work (in parallel) to assemble `K` and `r`
work!(assembler, buffer);

# As a side-note; the `StationaryFourier` element doesn't include a body load/source term contribution,
# such as the element in the heat flow tutorial. To add this contribution to the residual,
# we can use the `LoadHandler` with a `BodyLoad`:
lh = LoadHandler(dh)
add!(lh, BodyLoad(:u, 2, Returns(-1.0))) # rᵢ -= ∫ δuᵢ*1.0*dV 
apply!(r, lh, 0.0) # t = 0
