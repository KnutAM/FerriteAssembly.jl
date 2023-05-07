# # Plasticity
# This example shows how any material following the 
# [`MaterialModelsBase.jl`](https://github.com/KnutAM/MaterialModelsBase.jl)
# interface can be assembled with `FerriteAssembly.jl`. 
# [`element_routine!`](@ref FerriteAssembly.element_routine!(Ke, re, state::Vector{<:FerriteAssembly.MMB.AbstractMaterialState}, ae, material::FerriteAssembly.MMB.AbstractMaterial, cellvalues::CellVectorValues, dh_fh, Δt, cb))
# has already been implemented for this material. This example shows how to 
# * Assemble materials with state variables
# * Update state variables for the next time step 
# 
# ## Implementation
# We start by the required packages
using Tensors, MaterialModelsBase, Ferrite, FerriteAssembly

# And then we define the `material_response` for the plasticity material in
# [`J2Plasticity.jl`](J2Plasticity.jl), which is basically the same as in 
# [`Ferrite.jl`'s plasticity example](https://ferrite-fem.github.io/Ferrite.jl/stable/examples/plasticity/)
include("J2Plasticity.jl");

# ## Assembly
# With all required functions defined, we can now setup and assemble the finite element problem 
material = J2Plasticity(200.0e9, 0.3, 200.0e6, 10.0e9);
grid = generate_grid(Tetrahedron, (20,2,4), zero(Vec{3}), Vec((10.0,1.0,1.0)));
cellvalues = CellVectorValues(
    QuadratureRule{3,RefTetrahedron}(2), Lagrange{3, RefTetrahedron, 1}());
dh = DofHandler(grid); add!(dh, :u, 3); close!(dh); # Create dofhandler
K = create_sparsity_pattern(dh);
r = zeros(ndofs(dh));

# Using the `setup_assembly` function, 
buffer, states_old, states_new = setup_assembly(dh, material, cellvalues)
# we setup the `buffer`, old state variables, and new state variables. 
# The state variables are created via the [`create_cell_state`](@ref) 
# function that is already defined for `MaterialModelsBase.AbstractMaterial`

# We can now just provide an initial guess for the degree of freedom vector,`a`,
# and do the assembly
a = zeros(ndofs(dh))
doassemble!(K, r, states_new, states_old, buffer; a=a);

# If we would have a rate-dependent material, such that the time increment mattered,
# we can also supply that (but that is not required in this example)
doassemble!(K, r, states_new, states_old, buffer; a=a, Δt=1.0);

# In a full FE-program we iterate until convergence to find `a`. When converged,
# we go to the next time step, and would like to set the old state equal to the 
# current state, which we can do by calling 
update_states!(states_old, states_new)

# If we would like to access the states in any cell, states_old and states_new can be 
# indexed with the cell number. 