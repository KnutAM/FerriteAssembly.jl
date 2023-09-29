using Ferrite, FerriteAssembly
grid = generate_grid(Quadrilateral, (20, 20))
ip = Lagrange{RefQuadrilateral,1}()
dh = DofHandler(grid); add!(dh, :u, ip); close!(dh)
cellvalues = CellValues(QuadratureRule{2, RefCube}(2), ip);
K = create_sparsity_pattern(dh)
r = zeros(ndofs(dh));

struct ThermalMaterial
    k::Float64 # Thermal conductivity
    f::Float64 # Volumetric heat source
end
material = ThermalMaterial(1.0, 1.0);

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

ch = ConstraintHandler(dh)
faceset = union((getfaceset(grid,k) for k in ("left", "right", "bottom", "top"))...)
add!(ch, Dirichlet(:u, faceset, Returns(0.0)))
close!(ch);
apply_zero!(K, r, ch)

a = -K\r
vtk_grid("heat_equation", grid) do vtk
    vtk_point_data(vtk, dh, a)
end;

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
