# # Plasticity
# This example shows how any material following the 
# [`MaterialModelsBase.jl`](https://github.com/KnutAM/MaterialModelsBase.jl)
# interface can be assembled with `FerriteAssembly.jl`. 
# We start by the required packages
using Tensors, MaterialModelsBase, Ferrite, FerriteAssembly

include("J2Plasticity.jl")

include("MaterialModelsBaseElement.jl")

# ## Assembly
# With all required functions defined, we can now setup and assemble the finite element problem 
material = J2Plasticity(200.0e9, 0.3, 200.0e6, 10.0e9);
grid = generate_grid(Tetrahedron, (20,2,4), zero(Vec{3}), Vec((10.0,1.0,1.0)));
cellvalues = CellVectorValues(
    QuadratureRule{3,RefTetrahedron}(2), Lagrange{3, RefTetrahedron, 1}());
dh = DofHandler(grid); push!(dh, :u, 3); close!(dh); # Create dofhandler
K = create_sparsity_pattern(dh);
r = zeros(ndofs(dh));

# Using the `create_states` function, we can easily create the storage of state variables 
# suitable for the chosen dof handler and the given material (via the `create_cell_state`
# function that we defined earlier)
states = create_states(dh, material, cellvalues);

# And then we create the buffer for saving cell-related variables
buffer = setup_cellbuffer(dh, cellvalues, material, nothing, get_cache(material));

# Finally, with an initial guess of displacements, `a`, 
# we can create the assembler and do the assembly
a = zeros(ndofs(dh))
assembler = start_assemble(K,r)
doassemble!(assembler, buffer, states, dh, a);