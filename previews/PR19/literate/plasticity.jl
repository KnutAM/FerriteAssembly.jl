# # Plasticity
# This example shows how to assemble a plastic material following the 
# [`MaterialModelsBase.jl`](https://github.com/KnutAM/MaterialModelsBase.jl)
# interface with `FerriteAssembly.jl`. The 
# [`element_routine!`](@ref FerriteAssembly.element_routine!(Ke, re, state::Vector{<:MMB.AbstractMaterialState}, ae, material::MMB.AbstractMaterial, cellvalues::CellVectorValues, buffer))
# function implementation for any `MaterialModelsBase.AbstractMaterial` is already defined. \\
# Specifically, how to
# * Assemble materials with state variables
# * Giving the time increment for rate dependent elements
# * Update state variables for the next time step 
# 
# ## Material modeling 
# We start by the required packages
using Tensors, MaterialModelsBase, Ferrite, FerriteAssembly

# And then we define the `material_response` for the plasticity material in
# [`J2Plasticity.jl`](J2Plasticity.jl), which is basically the same as in 
# [`Ferrite.jl`'s plasticity example](https://ferrite-fem.github.io/Ferrite.jl/stable/examples/plasticity/)
include("J2Plasticity.jl");

# ## Standard `Ferrite.jl` setup
# With all required functions defined, we can now setup and assemble the finite element problem 
material = J2Plasticity(200.0e9, 0.3, 200.0e6, 10.0e9);
grid = generate_grid(Tetrahedron, (20,2,4), zero(Vec{3}), Vec((10.0,1.0,1.0)));
cellvalues = CellVectorValues(
    QuadratureRule{3,RefTetrahedron}(2), Lagrange{3, RefTetrahedron, 1}());
dh = DofHandler(grid); add!(dh, :u, 3); close!(dh); # Create dofhandler
K = create_sparsity_pattern(dh);
r = zeros(ndofs(dh));

# ## Setting up the assembly
# Using the `setup_assembly` function, 
buffer = setup_domainbuffer(DomainSpec(dh, material, cellvalues));
# we setup the `buffer`, old state variables, and new state variables. 
# The state variables are created via the [`create_cell_state`](@ref FerriteAssembly.create_cell_state) 
# function that is already defined for `MaterialModelsBase.AbstractMaterial`

# ## Doing the assembly
# We can now just provide an initial guess for the degree of freedom vector,`a`,
# and do the assembly
a = zeros(ndofs(dh))
assembler = start_assemble(K, r)
work!(assembler, buffer; a=a);

# If we would have a rate-dependent material, such that the time increment mattered,
# we can also supply that (but that is not required in this example)
set_time_increment!(buffer, 1.0)
assembler = start_assemble(K, r)
work!(assembler, buffer; a=a);

# ## Updating state variables
# In a full FE-program we iterate until convergence to find `a`. When converged,
# we go to the next time step, and would like to set the old state equal to the 
# current state, which we can do by calling 
update_states!(buffer);

# If we would like to access the states in any cell, we can request that from the buffer
cell_state = FerriteAssembly.get_state(buffer, 1)
display(typeof(cell_state))