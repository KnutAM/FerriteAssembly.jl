# # Multiple materials
# *How to assemble with different materials on different parts of the grid*\
# Specifically, how to
# * Setup multiple domains
# * Run threaded assembly
# 
# ## Material modeling 
# In addition to `J2Plasticity`, we setup a portion of the domain to use the `ElasticMaterial`
# that is also defined in [`J2Plasticity.jl`](J2Plasticity.jl) file.
using Tensors, MaterialModelsBase, Ferrite, FerriteAssembly
include("J2Plasticity.jl");

# ## Standard `Ferrite.jl` setup
# We start by setting up the 
grid = generate_grid(Tetrahedron, (20,2,4), zero(Vec{3}), Vec((10.0,1.0,1.0)))
cellvalues = CellVectorValues(
    QuadratureRule{3,RefTetrahedron}(2), Lagrange{3, RefTetrahedron, 1}())
dh = DofHandler(grid); add!(dh, :u, 3); close!(dh) # Create dofhandler
K = create_sparsity_pattern(dh)
r = zeros(ndofs(dh))
a = zeros(ndofs(dh));

# ## Setting up domains
# In order to setup a simulation with multiple domains, 
# we create one `DomainSpec` for each domain. 
# We start by creating the **elastic domain**,
addcellset!(grid, "elastic", x -> x[1] <= 5.0+eps())
cellset_el = getcellset(grid, "elastic")
material_el = ElasticMaterial(E=200.0e9, ν=0.3)
domain_el = DomainSpec(dh, material_el, cellvalues; set=cellset_el);

# followed by the **plastic domain**
cellset_pl = setdiff(1:getncells(grid), cellset_el)
material_pl = J2Plasticity(200.0e9, 0.3, 200.0e6, 10.0e9)
domain_pl = DomainSpec(dh, material_pl, cellvalues; set=cellset_pl);

# And then we can set up the assembly, where `threading=true` makes it multithreaded. 
buffers = setup_domainbuffers(Dict("elastic"=>domain_el, "plastic"=>domain_pl); threading=true);
# For multiple domains, `buffers::Dict{String}` with same keys as the dictionary given to the 
# setup function

# ## Doing the assembly
assembler = start_assemble(K, r)
work!(assembler, buffers; a=a);

# ## Updating state variables
# Updating the state variables after convergence in the current time step works as for single domains,
update_states!(buffers);

# Although the material behaviors are different, #src
# there are no differences in the responses as the displacements are zero   #src
# Hence, we can verify that we get the same stiffness in both cases     #src
# Running a test to be sure #src
using Test #hide
K_ref = create_sparsity_pattern(dh) #hide
r_ref = zeros(ndofs(dh)) #hide
a_ref = zeros(ndofs(dh)) #hide
buffer = setup_domainbuffer(DomainSpec(dh, material_el, cellvalues)) #hide
assembler_ref = start_assemble(K_ref, r_ref) #hide
work!(assembler_ref, buffer; a=a_ref) #hide
@test K ≈ K_ref    #hide
nothing;           #hide