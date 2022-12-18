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

materials = Dict(
    "elastic" => ElasticMaterial(E=200.0e9, ν=0.3),
    "plastic" => J2Plasticity(200.0e9, 0.3, 200.0e6, 10.0e9));
grid = generate_grid(Tetrahedron, (20,2,4), zero(Vec{3}), Vec((10.0,1.0,1.0)));
cellvalues = CellVectorValues(
    QuadratureRule{3,RefTetrahedron}(2), Lagrange{3, RefTetrahedron, 1}());
dh = DofHandler(grid); push!(dh, :u, 3); close!(dh); # Create dofhandler
K = create_sparsity_pattern(dh);
r = zeros(ndofs(dh));

addcellset!(grid, "elastic", x -> x[1] <= 5.0+eps())
addcellset!(grid, "plastic", setdiff(1:getncells(grid), getcellset(grid,"elastic")));

statefuns = Dict(key=>_->initial_material_state(mat) for (key,mat) in materials)
states = create_states(dh, statefuns, cellvalues);

caches = Dict(key=>get_cache(mat) for (key,mat) in materials)
buffers = setup_cellbuffer(dh, cellvalues, materials, nothing, caches);

a = zeros(ndofs(dh))
assembler = start_assemble(K, r)
doassemble!(assembler, buffers, states, dh, a);

using Test #hide
K_ref = create_sparsity_pattern(dh); #hide
r_ref = zeros(ndofs(dh)); #hide
states_ref = create_states(dh, _->initial_material_state(materials["elastic"]), cellvalues); #hide
a_ref = zeros(ndofs(dh)) #hide
assembler_ref = start_assemble(K_ref,r_ref) #hide
doassemble!(assembler_ref, buffers["elastic"], states_ref, dh, a_ref); #hide
@test K ≈ K_ref; #hide

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

