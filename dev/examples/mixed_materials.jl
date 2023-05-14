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
cellset_el = getcellset(grid, "elastic")
material_el = ElasticMaterial(E=200.0e9, ν=0.3)
domain_el = AssemblyDomain("elast", dh, material_el, cellvalues; cellset=cellset_el);

cellset_pl = setdiff(1:getncells(grid), cellset_el)
material_pl = J2Plasticity(200.0e9, 0.3, 200.0e6, 10.0e9)
domain_pl = AssemblyDomain("plast", dh, material_pl, cellvalues; cellset=cellset_pl);

buffers, new_states, old_states = setup_assembly([domain_el, domain_pl]; threading=true);

assembler = start_assemble(K, r)
doassemble!(assembler, new_states, buffers; a=a, old_states=old_states);

update_states!(old_states, new_states);

using Test #hide
K_ref = create_sparsity_pattern(dh) #hide
r_ref = zeros(ndofs(dh)) #hide
a_ref = zeros(ndofs(dh)) #hide
buffer, new_states, old_states = setup_assembly(dh, material_el, cellvalues) #hide
assembler_ref = start_assemble(K_ref, r_ref)
doassemble!(assembler_ref, new_states, buffer; a=a_ref, old_states=old_states) #hide
@test K ≈ K_ref    #hide
nothing;           #hide

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

