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

buffers, old_states, new_states = setup_assembly([elastic_domain, plastic_domain]; threading=true);

assembler = start_assemble(K, r)
doassemble!(assembler, new_states, buffers; a=a, old_states=old_states);

using Test #hide
K_ref = create_sparsity_pattern(dh) #hide
r_ref = zeros(ndofs(dh)) #hide
a_ref = zeros(ndofs(dh)) #hide
buffer, old_states, new_states = setup_assembly(dh, elastic_material, cellvalues) #hide
assembler_ref = start_assemble(K_ref, r_ref)
doassemble!(assembler_ref, new_states, buffer; a=a_ref, old_states=old_states) #hide
@test K ≈ K_ref    #hide

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

