using Ferrite, FerriteAssembly, BenchmarkTools
dh = DofHandler(generate_grid(Quadrilateral, (100, 100))); add!(dh, :u, 1); close!(dh)
cellvalues = CellScalarValues(QuadratureRule{2, RefCube}(2), Lagrange{2, RefCube, 1}());
K = create_sparsity_pattern(dh)
r = zeros(ndofs(dh));

struct ThermalMaterial
    k::Float64 # Thermal conductivity
    f::Float64 # Volumetric heat source
end
material = ThermalMaterial(1.0, 1.0)

function FerriteAssembly.element_routine!(Ke, re, state, ae,
        material::ThermalMaterial, cellvalues, cellbuffer
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

grid_domain = DomainSpec(dh, material, cellvalues)
buffer = setup_domainbuffer(grid_domain);

assembler = start_assemble(K, r);

work!(assembler, buffer);
K1 = deepcopy(K); #hide

threaded_buffer = setup_domainbuffer(grid_domain; threading=true);

assembler = start_assemble(K, r)
work!(assembler, threaded_buffer);
K2 = deepcopy(K); #hide

struct ThermalMaterialAD
    k::Float64 # Thermal conductivity
    f::Float64 # Volumetric heat source
end
material_ad = ThermalMaterialAD(1.0, 1.0)

function FerriteAssembly.element_residual!(re, state, ae,
        material::ThermalMaterialAD, cellvalues, cellbuffer
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

grid_domain_ad = DomainSpec(dh, material_ad, cellvalues)
buffer_ad = setup_domainbuffer(grid_domain_ad);

a = zeros(ndofs(dh))
assembler = start_assemble(K, r)
work!(assembler, buffer_ad; a=a);
K3 = deepcopy(K); #hide

if get(ENV, "CI", "false") == "true"        #hide
@btime work!($assembler, $buffer; a=$a)
@btime work!($assembler, $buffer_ad; a=$a)
end                                         #hide

buffer_ad2 = setup_domainbuffer(grid_domain_ad; autodiffbuffer=true)
if get(ENV, "CI", "false") == "true"        #hide
@btime work!($assembler, $buffer_ad2; a=$a)
end                                         #hide
                                                            #hide
assembler = start_assemble(K, r)                            #hide
work!(assembler, buffer_ad2; a=a);                          #hide
K4 = deepcopy(K);                                           #hide
                                                            #hide
using Test                                      #hide
@test K1 ≈ K2 ≈ K3 ≈ K4                         #hide
nothing;                                        #hide

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

