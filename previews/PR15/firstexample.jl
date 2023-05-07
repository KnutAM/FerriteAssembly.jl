using Ferrite, FerriteAssembly, BenchmarkTools
dh = DofHandler(generate_grid(Quadrilateral, (100, 100))); add!(dh, :u, 1); close!(dh)
cellvalues = CellScalarValues(QuadratureRule{2, RefCube}(2), Lagrange{2, RefCube, 1}());

struct ThermalMaterial end;

function FerriteAssembly.element_routine!(
        Ke::AbstractMatrix, re::AbstractVector, state, ae::AbstractVector,
        material::ThermalMaterial, cellvalues::CellScalarValues, buffer::FerriteAssembly.CellBuffer
        )
    n_basefuncs = getnbasefunctions(cellvalues)
    # Loop over quadrature points
    for q_point in 1:getnquadpoints(cellvalues)
        dΩ = getdetJdV(cellvalues, q_point)
        for i in 1:n_basefuncs
            δN  = shape_value(cellvalues, q_point, i)
            ∇δN = shape_gradient(cellvalues, q_point, i)
            # Add body load contribution to re
            re[i] += -δN * dΩ
            # Loop over trial shape functions
            for j in 1:n_basefuncs
                ∇N = shape_gradient(cellvalues, q_point, j)
                # Add contribution to Ke
                Ke[i, j] += (∇δN ⋅ ∇N) * dΩ
            end
        end
    end
end;

buffer, old_states, new_states = setup_assembly(dh, ThermalMaterial(), cellvalues)

K = create_sparsity_pattern(dh)
r = zeros(ndofs(dh));

doassemble!(K, r, new_states, old_states, buffer)

colors = create_coloring(dh.grid);

buffer2, _, _ = setup_assembly(dh, ThermalMaterial(), cellvalues; colors=colors);

doassemble!(K, r, new_states, old_states, buffer2);

struct ThermalMaterialAD end

function FerriteAssembly.element_residual!(
        re::AbstractVector, state, ae::AbstractVector,
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
            re[i] += (∇δN ⋅ ∇u - δN) * dΩ
        end
    end
end;

buffer_ad, old_states_ad, new_states_ad = setup_assembly(dh, ThermalMaterialAD(), cellvalues);

a = zeros(ndofs(dh))
doassemble!(K,r, new_states_ad, old_states_ad, buffer_ad; a=a);

@btime doassemble!($K, $r, $new_states, $old_states, $buffer; a=$a)
@btime doassemble!($K, $r, $new_states_ad, $old_states_ad, $buffer_ad; a=$a)

buffer_ad2, _, _ = setup_assembly(dh, ThermalMaterialAD(), cellvalues; autodiffbuffer=true)
@btime doassemble!($K, $r, $new_states_ad, $old_states_ad, $buffer_ad2; a=$a)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

