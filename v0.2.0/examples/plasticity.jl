using Tensors, MaterialModelsBase, Ferrite, FerriteAssembly

include("J2Plasticity.jl")

include("MaterialModelsBaseElement.jl")

material = J2Plasticity(200.0e9, 0.3, 200.0e6, 10.0e9);
grid = generate_grid(Tetrahedron, (20,2,4), zero(Vec{3}), Vec((10.0,1.0,1.0)));
cellvalues = CellVectorValues(
    QuadratureRule{3,RefTetrahedron}(2), Lagrange{3, RefTetrahedron, 1}());
dh = DofHandler(grid); push!(dh, :u, 3); close!(dh); # Create dofhandler
K = create_sparsity_pattern(dh);
r = zeros(ndofs(dh));

states = create_states(dh, material, cellvalues);

buffer = setup_cellbuffer(dh, cellvalues, material, nothing, get_cache(material));

a = zeros(ndofs(dh))
assembler = start_assemble(K,r)
doassemble!(assembler, buffer, states, dh, a);

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

