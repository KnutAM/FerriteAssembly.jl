@testset "heatequation" begin
    # Modified example from Ferrite.jl
    function get_grid()
        grid = generate_grid(Quadrilateral, (200, 200))
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

    weak = EE.WeakForm((δu, ∇δu, u, ∇u, u_dot, ∇u_dot) -> 1.0*(∇δu ⋅ ∇u) - δu*1.0)
    materials = (same=ThermalMaterial(), ad=ThermalMaterialAD(), weak=weak, mixed=Dict("A"=>ThermalMaterial(), "B"=>ThermalMaterialAD()))
    
    function setup_assembly_test(dh, material, cv; autodiff_cb=false, threaded=false)
        BufferType = threaded ? FerriteAssembly.ThreadedDomainBuffer : FerriteAssembly.DomainBuffer
        if isa(material, Dict) && isa(dh, DofHandler)
            setA, setB = (getcellset(dh.grid, name) for name in ("A", "B"))
            ad1 = DomainSpec(dh, material["A"], cv; set=setA)
            ad2 = DomainSpec(dh, material["B"], cv; set=setB)
            buffer = setup_domainbuffers(Dict("A"=>ad1, "B"=>ad2); autodiffbuffer=autodiff_cb, threading=threaded)
            @test isa(buffer, Dict{String,<:BufferType})
            @test isa(FerriteAssembly.get_old_state(buffer, "A"), Dict{Int})
            @test isa(FerriteAssembly.get_old_state(buffer, "B"), Dict{Int})
            return buffer
        elseif isa(material, Dict) && isa(dh, MixedDofHandler)
            sdh1 = FerriteAssembly.SubDofHandler(dh, dh.fieldhandlers[1])
            sdh2 = FerriteAssembly.SubDofHandler(dh, dh.fieldhandlers[2])
            set1 = getcellset(sdh1)
            setA, setB = (getcellset(dh.grid, name) for name in ("A", "B"))

            ad1 = DomainSpec(sdh1, material["A"], cv; set=intersect(setA, set1)) # sdh1A
            ad2 = DomainSpec(sdh1, material["B"], cv; set=intersect(setB, set1)) # sdh1B
            # For ad3 and ad4; add the full set to check correct intersection with sdh2's cellset internally. 
            ad3 = DomainSpec(sdh2, material["A"], cv; set=setA) # sdh2A 
            ad4 = DomainSpec(sdh2, material["B"], cv; set=setB) # sdh2B
            buffer = setup_domainbuffers(Dict("sdh1A"=>ad1, "sdh1B"=>ad2, "sdh2A"=>ad3, "sdh2B"=>ad4); autodiffbuffer=autodiff_cb, threading=threaded)
            @test isa(buffer, Dict{String,<:BufferType})
            @test isa(FerriteAssembly.get_old_state(buffer, "sdh1A"), Dict{Int})
            return buffer
        elseif isa(dh, MixedDofHandler)
            sdh1 = FerriteAssembly.SubDofHandler(dh, dh.fieldhandlers[1])
            sdh2 = FerriteAssembly.SubDofHandler(dh, dh.fieldhandlers[2])
            set1 = getcellset(sdh1); set2 = getcellset(sdh2)
            ad1 = DomainSpec(sdh1, material, cv; set=set1)
            ad2 = DomainSpec(sdh2, material, cv; set=set2)
            buffer = setup_domainbuffers(Dict("sdh1"=>ad1, "sdh2"=>ad2); autodiffbuffer=autodiff_cb, threading=threaded)
            @test isa(buffer, Dict{String,<:BufferType})
            @test isa(FerriteAssembly.get_old_state(buffer, "sdh1"), Dict{Int})
            return buffer
        else
            buffer = setup_domainbuffer(DomainSpec(dh, material, cv); autodiffbuffer=autodiff_cb, threading=threaded)
            @test isa(buffer, BufferType)
            @test isa(FerriteAssembly.get_old_state(buffer), Dict{Int})
            return buffer
        end
    end
    
    for DH in (DofHandler, MixedDofHandler)
        for mattype in (:same, :ad, :mixed, :weak)
            material = materials[mattype]
            
            cv, K, dh = setup_heatequation(DH)
            for scaling in (FerriteAssembly.NoScaling(), ElementResidualScaling(dh, 1))
                r = zeros(ndofs(dh))
                a = mattype==:same ? nothing : copy(r)  # If AD, dofs required

                @testset "$DH, $mattype, sequential" begin
                    autdiff_cbs = isa(material,ThermalMaterial) ? (false,) : (false, true)
                    for autodiff_cb in autdiff_cbs
                        fill!(K, 0); 
                        r .= rand(length(r)) # To ensure that it is actually changed
                        reset_scaling!(scaling)
                        buffer = setup_assembly_test(dh, material, cv; autodiff_cb=autodiff_cb)
                        ferrite_assembler = start_assemble(K, r)
                        assembler = isa(scaling, FerriteAssembly.NoScaling) ? ferrite_assembler : FerriteAssembly.KeReAssembler(ferrite_assembler; scaling=scaling)
                        work!(assembler, buffer; a=a)
                        isa(scaling, ElementResidualScaling) && @test scaling.factors[:u] ≈ sum(abs, r)  # As we use the 1-norm and all r's have the same sign
                        @test K_ref ≈ K 
                        @test r_ref ≈ r
                        if mattype == :ad    
                            # Assemble only r (requires element_residual!)
                            r .= rand(length(r)) # To ensure that it is both reset and then changed during assembly
                            reset_scaling!(scaling)
                            assembler = FerriteAssembly.ReAssembler(r; scaling=scaling)
                            work!(assembler, buffer; a=a)
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
                        ferrite_assembler = start_assemble(K, r)
                        assembler = isa(scaling, FerriteAssembly.NoScaling) ? ferrite_assembler : FerriteAssembly.KeReAssembler(ferrite_assembler; scaling=scaling)
                        buffer = setup_assembly_test(dh, material, cv; autodiff_cb=autodiff_cb, threaded=true)
                        # Quick check that test script works and that it is actually colored
                        TDB = FerriteAssembly.ThreadedDomainBuffer
                        @test isa(buffer, Union{Dict{String,<:TDB}, TDB})

                        work!(assembler, buffer; a=a)
                        if isa(scaling, ElementResidualScaling)
                            @test scaling.factors[:u] ≈ sum(abs, r)
                        end
                        @test K_ref ≈ K
                        @test r_ref ≈ r
                        
                        if mattype == :ad
                            r .= rand(length(r)) # To ensure that it is actually changed
                            reset_scaling!(scaling)
                            assembler = FerriteAssembly.ReAssembler(r; scaling=scaling)
                            work!(assembler, buffer; a=a)
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
