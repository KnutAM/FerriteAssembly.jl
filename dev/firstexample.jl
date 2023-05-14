using Ferrite, FerriteAssembly, BenchmarkTools
dh = DofHandler(generate_grid(Quadrilateral, (100, 100))); add!(dh, :u, 1); close!(dh)
cellvalues = CellScalarValues(QuadratureRule{2, RefCube}(2), Lagrange{2, RefCube, 1}());

struct ThermalMaterial
    k::Float64 # Thermal conductivity
    f::Float64 # Volumetric heat source
end

function FerriteAssembly.element_routine!(Ke, re, state, ae,
        material::ThermalMaterial, cellvalues, buffer
        )
    n_basefuncs = getnbasefunctions(cellvalues)
    # Loop over quadrature points
    for q_point in 1:getnquadpoints(cellvalues)
        dΩ = getdetJdV(cellvalues, q_point)
        for i in 1:n_basefuncs
            δN  = shape_value(cellvalues, q_point, i)
            ∇δN = shape_gradient(cellvalues, q_point, i)
            # Add body load contribution to re
            re[i] += -material.f*δN * dΩ
            # Loop over trial shape functions
            for j in 1:n_basefuncs
                ∇N = shape_gradient(cellvalues, q_point, j)
                # Add contribution to Ke
                Ke[i, j] += material.k*(∇δN ⋅ ∇N) * dΩ
            end
        end
    end
end;

material = ThermalMaterial(1.0, 1.0)
buffer, states = setup_assembly(dh, material, cellvalues);

K = create_sparsity_pattern(dh)
r = zeros(ndofs(dh));

assembler = start_assemble(K, r)
doassemble!(assembler, states, buffer);
K1 = deepcopy(K); #hide

threaded_buffer, _ = setup_assembly(dh, ThermalMaterial(1.0, 1.0), cellvalues; threading=true);

assembler = start_assemble(K, r)
doassemble!(assembler, states, threaded_buffer);
K2 = deepcopy(K); #hide

struct ThermalMaterialAD
    k::Float64 # Thermal conductivity
    f::Float64 # Volumetric heat source
end

function FerriteAssembly.element_residual!(re, state, ae,
        material::ThermalMaterialAD, cellvalues, buffer
        )
    n_basefuncs = getnbasefunctions(cellvalues)
    # Loop over quadrature points
    for q_point in 1:getnquadpoints(cellvalues)
        # Get the quadrature weight
        dΩ = getdetJdV(cellvalues, q_point)
        ∇u = function_gradient(cellvalues, q_point, ae)
        # Loop over test shape functions
        for i in 1:n_basefuncs
            δN  = shape_value(cellvalues, q_point, i)
            ∇δN = shape_gradient(cellvalues, q_point, i)
            # re = fint - fext
            re[i] += (material.k*(∇δN ⋅ ∇u) - material.f*δN) * dΩ
        end
    end
end;

material_ad = ThermalMaterialAD(1.0, 1.0)
buffer_ad, states_ad = setup_assembly(dh, material_ad, cellvalues);

a = zeros(ndofs(dh))
assembler = start_assemble(K, r)
doassemble!(assembler, states_ad, buffer_ad; a=a);
K3 = deepcopy(K); #hide

@btime doassemble!($assembler, $states, $buffer; a=$a)
@btime doassemble!($assembler, $states_ad, $buffer_ad; a=$a)

buffer_ad2, _, _ = setup_assembly(dh, material_ad, cellvalues; autodiffbuffer=true)
@btime doassemble!($assembler, $states_ad, $buffer_ad2; a=$a)
                                                            #hide
assembler = start_assemble(K, r)                            #hide
doassemble!(assembler, states_ad, buffer_ad2; a=a);         #hide
K4 = deepcopy(K);                                           #hide
                                                            #hide
using Test                                      #hide
@test K1 ≈ K2 ≈ K3 ≈ K4                         #hide
nothing;                                        #hide

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

