using Ferrite, FerriteAssembly, BenchmarkTools
dh = DofHandler(generate_grid(Quadrilateral, (100, 100))); push!(dh, :u, 1); close!(dh)
cellvalues = CellScalarValues(QuadratureRule{2, RefCube}(2), Lagrange{2, RefCube, 1}());

struct ThermalMaterial end;

function FerriteAssembly.element_routine!(
    Ke::AbstractMatrix, re::AbstractVector, state,
    ae::AbstractVector, material::ThermalMaterial, cellvalues::CellScalarValues,
    dh_fh, Δt, buffer::CellBuffer
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

K = create_sparsity_pattern(dh)
r = zeros(ndofs(dh));

states = create_states(dh);

cellbuffer = CellBuffer(dh, cellvalues, ThermalMaterial());

assembler = start_assemble(K,r)
doassemble!(assembler, cellbuffer, states, dh);

colors = create_coloring(dh.grid);

cellbuffers = create_threaded_CellBuffers(CellBuffer(dh, cellvalues, ThermalMaterial()))
assemblers = create_threaded_assemblers(K, r);

doassemble!(assemblers, cellbuffers, states, dh, colors);

struct ThermalMaterialAD end

function FerriteAssembly.element_residual!(
    re::AbstractVector, state,
    ae::AbstractVector, material::ThermalMaterialAD, cellvalues,
    dh_fh::Union{DofHandler,FieldHandler}, Δt, buffer
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

a = zeros(ndofs(dh))
cellbuffer2 = CellBuffer(dh, cellvalues, ThermalMaterialAD())
assembler = start_assemble(K,r)
doassemble!(assembler, cellbuffer2, states, dh, a);

@btime doassemble!(assembler, $cellbuffer, $states, $dh) setup=(assembler=start_assemble(K,r));
@btime doassemble!(assembler, $cellbuffer2, $states, $dh) setup=(assembler=start_assemble(K,r));

cellbuffer3 = AutoDiffCellBuffer(states, dh, cellvalues, ThermalMaterialAD())
@btime doassemble!(assembler, $cellbuffer3, $states, $dh) setup=(assembler=start_assemble(K,r));

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

