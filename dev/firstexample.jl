using Ferrite, FerriteAssembly
dh = DofHandler(generate_grid(Quadrilateral, (20, 20))); push!(dh, :u, 1); close!(dh)
cellvalues = CellScalarValues(QuadratureRule{2, RefCube}(2), Lagrange{2, RefCube, 1}())

struct ThermalMaterial end

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
            δu  = shape_value(cellvalues, q_point, i)
            ∇δu = shape_gradient(cellvalues, q_point, i)
            # Add body load contribution to re
            re[i] += -δu * dΩ
            # Loop over trial shape functions
            for j in 1:n_basefuncs
                ∇u = shape_gradient(cellvalues, q_point, j)
                # Add contribution to Ke
                Ke[i, j] += (∇δu ⋅ ∇u) * dΩ
            end
        end
    end
end

K = create_sparsity_pattern(dh)
r = zeros(ndofs(dh))

states = create_states(dh)

cellbuffer = CellBuffer(dh, cellvalues, ThermalMaterial())

assembler = start_assemble(K,r)
doassemble!(assembler, cellbuffer, states, dh)

colors = create_coloring(dh.grid)

cellbuffers = create_threaded_CellBuffers(CellBuffer(dh, cellvalues, ThermalMaterial()))
assemblers = create_threaded_assemblers(K, r)

doassemble!(assemblers, cellbuffers, states, colors, dh)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

