using Tensors, MaterialModelsBase, Ferrite, FerriteAssembly
include("J2Plasticity.jl");

grid = generate_grid(Tetrahedron, (20,2,4), zero(Vec{3}), Vec((10.0,1.0,1.0)))
cellvalues = CellVectorValues(
    QuadratureRule{3,RefTetrahedron}(2), Lagrange{3, RefTetrahedron, 1}())
dh = DofHandler(grid); add!(dh, :u, 3); close!(dh) # Create dofhandler
K = create_sparsity_pattern(dh)
r = zeros(ndofs(dh))
a = zeros(ndofs(dh));

addcellset!(grid, "elastic", x -> x[1] <= 5.0+eps())
elastic_material = ElasticMaterial(E=200.0e9, ν=0.3)
elastic_domain = AssemblyDomain("elast", dh, elastic_material, cellvalues; cellset=getcellset(grid, "elastic"));

plastic_cellset = setdiff(1:getncells(grid), getcellset(grid,"elastic"))
plastic_material = J2Plasticity(200.0e9, 0.3, 200.0e6, 10.0e9)
plastic_domain = AssemblyDomain("plast", dh, plastic_material, cellvalues; cellset=plastic_cellset);

colors = create_coloring(grid);
buffers, states_old, states_new = setup_assembly([elastic_domain, plastic_domain]; colors=colors);

doassemble!(K, r, states_new, states_old, buffers, a=a);

using Test #hide
K_ref = create_sparsity_pattern(dh) #hide
r_ref = zeros(ndofs(dh)) #hide
a_ref = zeros(ndofs(dh)) #hide
buffer, states_old, states_new = setup_assembly(dh, elastic_material, cellvalues) #hide
doassemble!(K_ref, r_ref, states_old, states_new, buffer; a=a_ref) #hide
@test K ≈ K_ref;    #hide

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

