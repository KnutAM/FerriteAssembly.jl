# # Multiple materials
# This example shows how to use two different materials on a 
# grid with the same cells everywhere. The key with this example is how
# * Setup multiple [`AssemblyDomain`](@ref)s
# * Run threaded assembly
# 
# In addition to `J2Plasticity`, we setup a portion of the domain to use the `ElasticMaterial`
# that is also defined in [`J2Plasticity.jl`](J2Plasticity.jl) file.
using Tensors, MaterialModelsBase, Ferrite, FerriteAssembly
include("J2Plasticity.jl");

# We start by setting up the 
grid = generate_grid(Tetrahedron, (20,2,4), zero(Vec{3}), Vec((10.0,1.0,1.0)))
cellvalues = CellVectorValues(
    QuadratureRule{3,RefTetrahedron}(2), Lagrange{3, RefTetrahedron, 1}())
dh = DofHandler(grid); add!(dh, :u, 3); close!(dh) # Create dofhandler
K = create_sparsity_pattern(dh)
r = zeros(ndofs(dh))
a = zeros(ndofs(dh));

# In order to setup a simulation with multiple domains, we must use the 
# [`AssemblyDomain`](@ref) structure to setup the simulation.
# We start by creating the elastic domain,
addcellset!(grid, "elastic", x -> x[1] <= 5.0+eps())
elastic_material = ElasticMaterial(E=200.0e9, ν=0.3)
elastic_domain = AssemblyDomain("elast", dh, elastic_material, cellvalues; cellset=getcellset(grid, "elastic"));

# And then create the plastic domain
plastic_cellset = setdiff(1:getncells(grid), getcellset(grid,"elastic"))
plastic_material = J2Plasticity(200.0e9, 0.3, 200.0e6, 10.0e9)
plastic_domain = AssemblyDomain("plast", dh, plastic_material, cellvalues; cellset=plastic_cellset);

# We can then setup the assembly, and in this case we would like to do the assembly threaded,
# so we need to color the domain.
colors = create_coloring(grid);
buffers, states_old, states_new = setup_assembly([elastic_domain, plastic_domain]; colors=colors);
# In this case, `buffers`, `states_old`, and `states_new` are `Dict{String}` with keys according to the 
# names given to each `AssemblyDomain`. This is important for postprocessing, but for doing assembly,
# these can be passed directly:
doassemble!(K, r, states_new, states_old, buffers, a=a);

# Although the material behaviors are different, #src
# there are no differences in the responses as the displacements are zero   #src
# Hence, we can verify that we get the same stiffness in both cases     #src
# Running a test to be sure #src
using Test #hide
K_ref = create_sparsity_pattern(dh) #hide
r_ref = zeros(ndofs(dh)) #hide
a_ref = zeros(ndofs(dh)) #hide
buffer, states_old, states_new = setup_assembly(dh, elastic_material, cellvalues) #hide
doassemble!(K_ref, r_ref, states_old, states_new, buffer; a=a_ref) #hide
@test K ≈ K_ref;    #hide
