# Modified example from Ferrite.jl
function setup_heatequation()
    grid = generate_grid(Quadrilateral, (20, 20));
    dim = 2
    ip = Lagrange{dim, RefCube, 1}()
    qr = QuadratureRule{dim, RefCube}(2)
    cellvalues = CellScalarValues(qr, ip);

    dh = DofHandler(grid)
    push!(dh, :u, 1)
    close!(dh);

    ch = ConstraintHandler(dh);
    
    ∂Ω = union(
        getfaceset(grid, "left"),
        getfaceset(grid, "right"),
        getfaceset(grid, "top"),
        getfaceset(grid, "bottom"),
    );

    dbc = Dirichlet(:u, ∂Ω, (x, t) -> 0)
    add!(ch, dbc);

    close!(ch)
    update!(ch, 0.0);

    K = create_sparsity_pattern(dh)

    return cellvalues, K, dh
end

function assemble_element!(Ke::Matrix, fe::Vector, cellvalues::CellScalarValues)
    n_basefuncs = getnbasefunctions(cellvalues)
    # Reset to 0
    fill!(Ke, 0)
    fill!(fe, 0)
    # Loop over quadrature points
    for q_point in 1:getnquadpoints(cellvalues)
        # Get the quadrature weight
        dΩ = getdetJdV(cellvalues, q_point)
        # Loop over test shape functions
        for i in 1:n_basefuncs
            δu  = shape_value(cellvalues, q_point, i)
            ∇δu = shape_gradient(cellvalues, q_point, i)
            # Add contribution to fe
            fe[i] += δu * dΩ
            # Loop over trial shape functions
            for j in 1:n_basefuncs
                ∇u = shape_gradient(cellvalues, q_point, j)
                # Add contribution to Ke
                Ke[i, j] += (∇δu ⋅ ∇u) * dΩ
            end
        end
    end
    return Ke, fe
end

function assemble_global(cellvalues::CellScalarValues, K::SparseMatrixCSC, dh::DofHandler)
    # Allocate the element stiffness matrix and element force vector
    n_basefuncs = getnbasefunctions(cellvalues)
    Ke = zeros(n_basefuncs, n_basefuncs)
    fe = zeros(n_basefuncs)
    # Allocate global force vector f
    f = zeros(ndofs(dh))
    # Create an assembler
    assembler = start_assemble(K, f)
    # Loop over all cels
    for cell in CellIterator(dh)
        # Reinitialize cellvalues for this cell
        reinit!(cellvalues, cell)
        # Compute element contribution
        assemble_element!(Ke, fe, cellvalues)
        # Assemble Ke and fe into K and f
        assemble!(assembler, celldofs(cell), Ke, fe)
    end
    return K, f
end

struct ThermalMaterial end  # For dispatch only, material parameters hard coded in original example

function FerriteAssembly.element_routine!(
    Ke::AbstractMatrix, re::AbstractVector, 
    ae_new::AbstractVector, ae_old::AbstractVector,
    state::AbstractVector, material::ThermalMaterial, 
    cellvalues::CellScalarValues, 
    dh_fh::Union{DofHandler,FieldHandler}, Δt, materialcache
    )
    assemble_element!(Ke, re, cellvalues)
    re .*= -1   # re = fint-fext
    return nothing
end

@testset "heatequation" begin
    # Following Ferrite.jl's example:
    cv_ref, K_ref, dh_ref = setup_heatequation()
    K_ref, f_ref = assemble_global(cv_ref, K_ref, dh_ref);

    # Using FerriteAssembly:
    cv, K, dh = setup_heatequation()
    
    states = [[nothing for _ in 1:getnquadpoints(cv)] for _ in 1:getncells(dh.grid)]
    material = ThermalMaterial()
    a=zeros(ndofs(dh)); aold=copy(a);
    Δt=1.0
    r = zeros(ndofs(dh))
    cache = CellCache(dh)
    doassemble!(K, r, a, aold, states, dh, cv, material, Δt, cache)

    @test K_ref ≈ K 
    @test f_ref ≈ -r
end