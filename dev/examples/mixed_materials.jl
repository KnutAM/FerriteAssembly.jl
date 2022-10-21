using Tensors, MaterialModelsBase, Ferrite, FerriteAssembly
include("J2Plasticity.jl");
include("MaterialModelsBaseElement.jl");

struct ElasticMaterial{T<:SymmetricTensor{4,3}} <: AbstractMaterial
    D::T
end
function ElasticMaterial(;E, ν)
    D, _, _ = elastic_conversion(E, ν)
    return ElasticMaterial(D)
end;

function MaterialModelsBase.material_response(
    material::ElasticMaterial, ϵ::SymmetricTensor{2,3}, state,
    Δt, cache=get_cache(material), args...; kwargs...)
    σ = material.D ⊡ ϵ
    return σ, material.D, NoMaterialState()
end;

material_pl = J2Plasticity(200.0e9, 0.3, 200.0e6, 10.0e9);
material_el = ElasticMaterial(E=200.0e9, ν=0.3);
grid = generate_grid(Tetrahedron, (20,2,4), zero(Vec{3}), Vec((10.0,1.0,1.0)));
cellvalues = CellVectorValues(
    QuadratureRule{3,RefTetrahedron}(2), Lagrange{3, RefTetrahedron, 1}());
dh = DofHandler(grid); push!(dh, :u, 3); close!(dh); # Create dofhandler
K = create_sparsity_pattern(dh);
r = zeros(ndofs(dh));

addcellset!(grid, "elastic", x -> x[1] <= 5.0+eps())
set_el = getcellset(grid, "elastic");
addcellset!(grid, "plastic", setdiff(1:getncells(grid), set_el));
set_pl = getcellset(grid, "plastic");

states_el = create_states(dh, _->initial_material_state(material_el), cellvalues, set_el);
states_pl = create_states(dh, _->initial_material_state(material_pl), cellvalues, set_pl);

buffer_el = CellBuffer(dh, cellvalues, material_el, nothing, get_cache(material_el));
buffer_pl = CellBuffer(dh, cellvalues, material_pl, nothing, get_cache(material_pl));

a = zeros(ndofs(dh))
assembler = start_assemble(K, r)
doassemble!(assembler, buffer_el, states_el, dh, a; cellset=set_el);
doassemble!(assembler, buffer_pl, states_pl, dh, a; cellset=set_pl);

using Test #hide
K_ref = create_sparsity_pattern(dh); #hide
r_ref = zeros(ndofs(dh)); #hide
states = create_states(dh, _->initial_material_state(material_el), cellvalues); #hide
a_ref = zeros(ndofs(dh)) #hide
assembler_ref = start_assemble(K_ref,r_ref) #hide
doassemble!(assembler_ref, buffer_el, states, dh, a_ref); #hide
@test K ≈ K_ref; #hide

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

