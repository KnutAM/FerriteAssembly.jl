# Modified example from Ferrite.jl
function get_grid()
    grid = generate_grid(Quadrilateral, (20, 20))
    # Add cellsets for testing different materials on the grid
    addcellset!(grid, "A", x->x[1]<0)
    addcellset!(grid, "B", setdiff(1:getncells(grid), getcellset(grid,"A")))
    return grid
end

function get_dh()
    grid = get_grid();
    dh = DofHandler(grid)
    push!(dh, :u, 1)
    close!(dh);
    return dh
end

function get_mdh(ip)
    grid = get_grid();
    dh = MixedDofHandler(grid)
    push!(dh, FieldHandler([Field(:u, ip, 1)], Set(collect(1:getncells(grid)))))
    close!(dh);
    return dh
end

function setup_heatequation(DH=DofHandler)
    dim = 2
    ip = Lagrange{dim, RefCube, 1}()
    qr = QuadratureRule{dim, RefCube}(2)
    cellvalues = CellScalarValues(qr, ip);
    dh = DH==MixedDofHandler ? get_mdh(ip) : get_dh()
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
    Ke::AbstractMatrix, re::AbstractVector, state,
    ae::AbstractVector, material::ThermalMaterial, cellvalues, 
    dh_fh::Union{DofHandler,FieldHandler}, Δt, buffer::CellBuffer
    )
    assemble_element!(Ke, re, cellvalues)
    re .*= -1   # re = fint-fext
    return nothing
end

@testset "heatequation" begin
    # Following Ferrite.jl's example get get reference solution
    cv_ref, K_ref, dh_ref = setup_heatequation()
    K_ref, f_ref = assemble_global(cv_ref, K_ref, dh_ref);
    r_ref = -f_ref

    materials = (same=ThermalMaterial(), mixed=Dict("A"=>ThermalMaterial(), "B"=>ThermalMaterial()))

    for DH in (DofHandler, MixedDofHandler)
        for mattype in (:same, :mixed)
            material = materials[mattype]

            cv, K, dh = setup_heatequation(DH)
            r = zeros(ndofs(dh))
            if mattype==:same
                states = create_states(dh)
            else
                states = create_states(dh, Dict(key=>Returns(nothing) for key in keys(material)))
            end

            @testset "$DH, $mattype, sequential" begin
                cellbuffer = CellBuffer(dh, cv, material)
                assembler = start_assemble(K, r)
                doassemble!(assembler, cellbuffer, states, dh)

                @test K_ref ≈ K 
                @test r_ref ≈ r
            end

            @testset "$DH, $mattype, threaded" begin
                cellbuffers = create_threaded_CellBuffers(CellBuffer(dh, cv, material))
                assemblers = create_threaded_assemblers(K, r)
                colors = create_coloring(dh.grid)
                
                doassemble!(assemblers, cellbuffers, states, dh, colors)
        
                @test K_ref ≈ K 
                @test r_ref ≈ r
                
            end
        end
    end

end