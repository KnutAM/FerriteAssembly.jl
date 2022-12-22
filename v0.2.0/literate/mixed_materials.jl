# # Multiple materials
# This example shows how to use two different materials on a 
# grid with the same cells everywhere. Hence, the only difference 
# is the type of material and the state variables. 
# We start by including the necessary packages, as well as the 
# [`J2Plasticity.jl`](J2Plasticity.jl) and 
# [`MaterialModelsBaseElement.jl`](MaterialModelsBaseElement.jl)
# files from the previous example. 
using Tensors, MaterialModelsBase, Ferrite, FerriteAssembly
include("J2Plasticity.jl");
include("MaterialModelsBaseElement.jl");

# Then, we also define an elastic material 
struct ElasticMaterial{T<:SymmetricTensor{4,3}} <: AbstractMaterial
    D::T
end
function ElasticMaterial(;E, ν)
    D, _, _ = elastic_conversion(E, ν)
    return ElasticMaterial(D)
end;

# This material requires not state, so `MaterialModelsBase.jl`'s `NoMaterialState`
# will be created by default. Hence, we only need to define the `material_response`
function MaterialModelsBase.material_response(
    material::ElasticMaterial, ϵ::SymmetricTensor{2,3}, state, 
    Δt, cache=get_cache(material), args...; kwargs...)
    σ = material.D ⊡ ϵ
    return σ, material.D, NoMaterialState()
end;

# The same `element_routine!` as before can be used, defined in 
# [`MaterialModelsBaseElement.jl`](MaterialModelsBaseElement.jl).
# Hence, we just need to define the materials, grid, cellvalues etc.
materials = Dict(
    "elastic" => ElasticMaterial(E=200.0e9, ν=0.3),
    "plastic" => J2Plasticity(200.0e9, 0.3, 200.0e6, 10.0e9));
grid = generate_grid(Tetrahedron, (20,2,4), zero(Vec{3}), Vec((10.0,1.0,1.0)));
cellvalues = CellVectorValues(
    QuadratureRule{3,RefTetrahedron}(2), Lagrange{3, RefTetrahedron, 1}());
dh = DofHandler(grid); push!(dh, :u, 3); close!(dh); # Create dofhandler
K = create_sparsity_pattern(dh);
r = zeros(ndofs(dh));

# We then define the cellsets where with the same key as each material in `materials`
addcellset!(grid, "elastic", x -> x[1] <= 5.0+eps())
addcellset!(grid, "plastic", setdiff(1:getncells(grid), getcellset(grid,"elastic")));

# The `create_states` function will then create the correct datastructure for the 
# states, (one Dict{Int} for each `material`)
states = create_states(dh, materials, cellvalues);

# If desired, we can create a cache based on `MaterialModelsBase`, even though not used in this 
# example (we could also skip passing this to `setup_cellbuffer`)
caches = Dict(key=>get_cache(mat) for (key,mat) in materials);

# Then we create the cell buffers
buffers = setup_cellbuffer(dh, cellvalues, materials, nothing, caches);

# And then we define the displacements and the assembler, before assembling
a = zeros(ndofs(dh))
assembler = start_assemble(K, r)
doassemble!(assembler, buffers, states, dh, a);

# Although the material behaviors are different, #src
# there are no differences in the responses as the displacements are zero   #src
# Hence, we can verify that we get the same stiffness in both cases     #src
# Running a test to be sure #src
using Test #hide
K_ref = create_sparsity_pattern(dh); #hide
r_ref = zeros(ndofs(dh)); #hide
states_ref = create_states(dh, materials["elastic"], cellvalues); #hide
a_ref = zeros(ndofs(dh)) #hide
assembler_ref = start_assemble(K_ref,r_ref) #hide
doassemble!(assembler_ref, buffers["elastic"], states_ref, dh, a_ref); #hide
@test K ≈ K_ref; #hide
