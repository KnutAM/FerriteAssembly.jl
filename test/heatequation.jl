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
        add!(dh, :u, 1)
        close!(dh);
        return dh
    end

    function get_mdh(ip)
        grid = get_grid();
        dh = MixedDofHandler(grid)
        n_half = getncells(grid)÷2
        add!(dh, FieldHandler([Field(:u, ip, 1)], Set(collect(1:n_half))))
        add!(dh, FieldHandler([Field(:u, ip, 1)], Set(collect((n_half+1):getncells(grid)))))
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
            Ke::AbstractMatrix, re::AbstractVector, state, ae::AbstractVector, 
            material::ThermalMaterial, cellvalues, buffer)
        assemble_element!(Ke, re, cellvalues)
        re .*= -1   # re = fint-fext
        return nothing
    end

    struct ThermalMaterialAD end

    function FerriteAssembly.element_residual!(
            re::AbstractVector, state, ae::AbstractVector, 
            material::ThermalMaterialAD, cellvalues, buffer)
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
    cellbuffer = FerriteAssembly.setup_cellbuffer(false, FerriteAssembly.SubDofHandler(dh), cv, mtrl, [nothing], (u=dof_range(dh, :u),), nothing, nothing)
    cellbuffer_ad = FerriteAssembly.AutoDiffCellBuffer(cellbuffer)
    ae = FerriteAssembly.get_ae(cellbuffer)
    re = FerriteAssembly.get_re(cellbuffer)
    Ke = FerriteAssembly.get_Ke(cellbuffer)
    assemble_element!(Ke, re, cv)
    Ke_ref, fe_ref = copy.((Ke, re))
    for cb in (cellbuffer, cellbuffer_ad)
        fill!.((ae,re,Ke), 0)
        FerriteAssembly.element_routine!(Ke, re, nothing, ae, mtrl, cv, cellbuffer)
        @test Ke ≈ Ke_ref 
        @test re ≈ -fe_ref # as ae=0
    end
    weak = EE.WeakForm((δu, ∇δu, u, ∇u, u_dot, ∇u_dot) -> 1.0*(∇δu ⋅ ∇u) - δu*1.0)
    materials = (same=ThermalMaterial(), ad=ThermalMaterialAD(), weak=weak, mixed=Dict("A"=>ThermalMaterial(), "B"=>ThermalMaterialAD()))
    
    function setup_assembly_test(dh, material, cv, scaling, autodiff_cb, colors=nothing)
        if isa(material, Dict) && isa(dh, DofHandler)
            setA, setB = (getcellset(dh.grid, name) for name in ("A", "B"))
            ad1 = FerriteAssembly.AssemblyDomain("A", dh, material["A"], cv; cellset=setA)
            ad2 = FerriteAssembly.AssemblyDomain("B", dh, material["B"], cv; cellset=setB)
            buffer, states_old, states_new = setup_assembly([ad1, ad2]; autodiffbuffer=autodiff_cb, scaling=scaling, colors=colors)
            @test isa(buffer, Dict{String,<:FerriteAssembly.DomainBuffer})
            @test isa(states_old, Dict{String,<:Dict{Int}})
            return buffer, states_old, states_new
        elseif isa(material, Dict) && isa(dh, MixedDofHandler)
            sdh1 = FerriteAssembly.SubDofHandler(dh, dh.fieldhandlers[1])
            sdh2 = FerriteAssembly.SubDofHandler(dh, dh.fieldhandlers[2])
            set1 = getcellset(sdh1)
            setA, setB = (getcellset(dh.grid, name) for name in ("A", "B"))

            ad1 = FerriteAssembly.AssemblyDomain("sdh1A", sdh1, material["A"], cv; cellset=intersect(setA, set1))
            ad2 = FerriteAssembly.AssemblyDomain("sdh1B", sdh1, material["B"], cv; cellset=intersect(setB, set1))
            # For ad3 and ad4; add the full set to check correct intersection with sdh2's cellset internally. 
            ad3 = FerriteAssembly.AssemblyDomain("sdh2A", sdh2, material["A"], cv; cellset=setA)
            ad4 = FerriteAssembly.AssemblyDomain("sdh2B", sdh2, material["B"], cv; cellset=setB)
            buffer, states_old, states_new = setup_assembly([ad1, ad2, ad3, ad4]; autodiffbuffer=autodiff_cb, scaling=scaling, colors=colors)
            @test isa(buffer, Dict{String,<:FerriteAssembly.DomainBuffer})
            @test isa(states_old, Dict{String,<:Dict{Int}})
            return buffer, states_old, states_new
        elseif isa(dh, MixedDofHandler)
            sdh1 = FerriteAssembly.SubDofHandler(dh, dh.fieldhandlers[1])
            sdh2 = FerriteAssembly.SubDofHandler(dh, dh.fieldhandlers[2])
            set1 = getcellset(sdh1); set2 = getcellset(sdh2)
            ad1 = FerriteAssembly.AssemblyDomain("sdh1", sdh1, material, cv; cellset=set1)
            ad2 = FerriteAssembly.AssemblyDomain("sdh2", sdh2, material, cv; cellset=set2)
            buffer, states_old, states_new = setup_assembly([ad1, ad2]; autodiffbuffer=autodiff_cb, scaling=scaling, colors=colors)
            @test isa(buffer, Dict{String,<:FerriteAssembly.DomainBuffer})
            @test isa(states_old, Dict{String,<:Dict{Int}})
            return buffer, states_old, states_new
        else
            buffer, states_old, states_new = setup_assembly(dh, material, cv; scaling=scaling, autodiffbuffer=autodiff_cb, colors=colors)
            @test isa(buffer, FerriteAssembly.DomainBuffer)
            @test isa(states_old, Dict{Int})
            return buffer, states_old, states_new
        end
    end
    
    for DH in (DofHandler, MixedDofHandler)
        for mattype in (:same, :ad, :mixed, :weak)
            material = materials[mattype]
            
            cv, K, dh = setup_heatequation(DH)
            for scaling in (nothing, ElementResidualScaling(dh, 1))
                r = zeros(ndofs(dh))
                a = mattype==:same ? nothing : copy(r)  # If AD, dofs required
                
                @testset "$DH, $mattype, sequential" begin
                    autdiff_cbs = isa(material,ThermalMaterial) ? (false,) : (false, true)
                    for autodiff_cb in autdiff_cbs
                        fill!(K, 0); 
                        r .= rand(length(r)) # To ensure that it is actually changed
                        reset_scaling!(scaling)
                        buffer, states_old, states_new = setup_assembly_test(dh, material, cv, scaling, autodiff_cb)
                        doassemble!(K, r, states_new, states_old, buffer; a=a)
                        isa(scaling, ElementResidualScaling) && @test scaling.factors[:u] ≈ sum(abs, r)  # As we use the 1-norm and all r's have the same sign
                        @test K_ref ≈ K 
                        @test r_ref ≈ r
                        if mattype == :ad    
                            # Assemble only r (requires element_residual!)
                            r .= rand(length(r)) # To ensure that it is actually changed
                            reset_scaling!(scaling)
                            doassemble!(r, states_new, states_old, buffer; a=a)
                            isa(scaling, ElementResidualScaling) && @test scaling.factors[:u] ≈ sum(abs, r)  # As we use the 1-norm and all r's have the same sign 
                            @test r_ref ≈ r
                        end
                    end
                end
                @testset "$DH, $mattype, threaded" begin
                    autdiff_cbs = isa(material,ThermalMaterial) ? (false,) : (false, true)
                    for autodiff_cb in autdiff_cbs
                        fill!(K, 0); 
                        r .= rand(length(r)) # To ensure that it is actually changed
                        reset_scaling!(scaling)
                        colors = create_coloring(dh.grid)
                        buffer, states_old, states_new = setup_assembly_test(dh, material, cv, scaling, autodiff_cb, colors)
                        # Quick check that test script works and that it is actually colored
                        if isa(buffer, Dict)
                            @test FerriteAssembly.iscolored(buffer)
                        else
                            @test isa(buffer, FerriteAssembly.DomainBuffer{true})
                        end
                        doassemble!(K, r, states_new, states_old, buffer; a=a)
                        if isa(scaling, ElementResidualScaling)
                            @test scaling.factors[:u] ≈ sum(abs, r)
                        end
                        @test K_ref ≈ K
                        @test r_ref ≈ r
                        
                        if mattype == :ad    
                            r .= rand(length(r)) # To ensure that it is actually changed
                            reset_scaling!(scaling)
                            doassemble!(r, states_new, states_old, buffer; a=a)
                            if isa(scaling, ElementResidualScaling)
                                @test scaling.factors[:u] ≈ sum(abs, r)
                            end
                            @test r_ref ≈ r
                        end
                    end
                end
            end
        end
    end
end
