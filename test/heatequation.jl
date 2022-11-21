@testset "heatequation" begin
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
        dh_fh::Union{DofHandler,FieldHandler}, Δt, buffer
        )
        assemble_element!(Ke, re, cellvalues)
        re .*= -1   # re = fint-fext
        return nothing
    end

    struct ThermalMaterialAD end

    function FerriteAssembly.element_residual!(
        re::AbstractVector, state,
        ae::AbstractVector, material::ThermalMaterialAD, cellvalues, 
        dh_fh::Union{DofHandler,FieldHandler}, Δt, buffer::CellBuffer
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
        return nothing
    end
    
    # Following Ferrite.jl's example get get reference solution
    cv_ref, K_ref, dh_ref = setup_heatequation()
    K_ref, f_ref = assemble_global(cv_ref, K_ref, dh_ref);
    r_ref = -f_ref

    # Make quick check that the AD material gives the same result directly without assembly 
    cv, _, dh = setup_heatequation(DofHandler)
    reinit!(cv, getcoordinates(dh.grid,1))
    mtrl = ThermalMaterialAD()
    cellbuffer = CellBuffer(dh, cv, mtrl)
    cellbuffer_ad = FerriteAssembly.AutoDiffCellBuffer([nothing,],dh,cv,mtrl)
    ae = FerriteAssembly.get_ae(cellbuffer)
    re = FerriteAssembly.get_re(cellbuffer)
    Ke = FerriteAssembly.get_Ke(cellbuffer)
    assemble_element!(Ke, re, cv)
    Ke_ref, fe_ref = copy.((Ke, re))
    for cb in (cellbuffer, cellbuffer_ad)
        fill!.((ae,re,Ke), 0)
        FerriteAssembly.element_routine!(Ke, re, nothing, ae, mtrl, cv, dh, 0.0, cellbuffer)
        @test Ke ≈ Ke_ref 
        @test re ≈ -fe_ref # as ae=0
    end
    materials = (same=ThermalMaterial(), ad=ThermalMaterialAD(), mixed=Dict("A"=>ThermalMaterial(), "B"=>ThermalMaterialAD()))

    for DH in (DofHandler, MixedDofHandler)
        for mattype in (:same, :ad, :mixed)
            material = materials[mattype]
            
            cv, K, dh = setup_heatequation(DH)
            r = zeros(ndofs(dh))
            a = mattype==:same ? nothing : copy(r)
            if isa(material, Dict)
                states = create_states(dh, Dict(key=>Returns(nothing) for key in keys(material)))
            else
                states = create_states(dh)
            end

            @testset "$DH, $mattype, sequential" begin
                cellbuffer = CellBuffer(dh, cv, material)
                cbs = isa(material,ThermalMaterial) ? (cellbuffer,) : (cellbuffer,FerriteAssembly.AutoDiffCellBuffer(states,dh,cv,material))
                for cb in cbs    
                    assembler = start_assemble(K, r)
                    doassemble!(assembler, cb, states, dh, a)

                    @test K_ref ≈ K 
                    @test r_ref ≈ r
                end
            end

            @testset "$DH, $mattype, threaded" begin
                cellbuffer = CellBuffer(dh, cv, material)
                cbs = isa(material,ThermalMaterial) ? (cellbuffer,) : (cellbuffer,FerriteAssembly.AutoDiffCellBuffer(states,dh,cv,material))
                for cb in cbs
                    cellbuffers = create_threaded_CellBuffers(cb)
                    assemblers = create_threaded_assemblers(K, r)
                    colors = create_coloring(dh.grid)
                    
                    doassemble!(assemblers, cellbuffers, states, dh, colors, a)
            
                    @test K_ref ≈ K 
                    @test r_ref ≈ r
                end
            end
        end
    end

end