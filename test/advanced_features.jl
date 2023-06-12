@testset "modify_material" begin
    function run_assembly(buffer, states, dh)
        K = create_sparsity_pattern(dh)
        r = zeros(ndofs(dh))
        assembler = start_assemble(K, r)
        doassemble!(assembler, states, buffer)
        return K
    end

    function test_modified_material(dh, b1, b2, s, k1, k2)
        K1 = run_assembly(b1, s, dh)
        K2 = run_assembly(b2, s, dh)
        @test K2*k1 ≈ K1*k2
        FerriteAssembly.modify_material!(m->FerriteAssembly.ExampleElements.StationaryFourier(k2), b1)
        K1mod = run_assembly(b1, s1, dh)
        @test K2 ≈ K1mod
        return K1
    end

    k1, k2 = (1.0, 2.0)
    grid = generate_grid(Hexahedron, (2,2,2))
    dh = DofHandler(grid); add!(dh, :u, 1); close!(dh)
    m1 = FerriteAssembly.ExampleElements.StationaryFourier(k1)
    m2 = FerriteAssembly.ExampleElements.StationaryFourier(k2)
    ip = Ferrite.default_interpolation(Hexahedron)
    cv = CellScalarValues(QuadratureRule{3,RefCube}(2), ip)
    # Sequential, single domain
    b1, s1, _ = setup_assembly(dh, m1, cv)
    b2, s2, _ = setup_assembly(dh, m2, cv)
    K_test1 = test_modified_material(dh, b1, b2, s1, k1, k2)
    # Threaded, single domain
    b1, s1, _ = setup_assembly(dh, m1, cv; threading=true)
    b2, s2, _ = setup_assembly(dh, m2, cv; threading=true)
    K_test2 = test_modified_material(dh, b1, b2, s1, k1, k2)
    # Threaded, multiple domains
    addcellset!(grid, "A", x -> x[1] > -eps())
    addcellset!(grid, "B", setdiff!(Set(1:getncells(grid)), getcellset(grid, "A")))
    ad1_A = AssemblyDomain("A", dh, m1, cv; cellset=getcellset(grid, "A"))
    ad1_B = AssemblyDomain("B", dh, m1, cv; cellset=getcellset(grid, "B"))
    ad2_A = AssemblyDomain("A", dh, m2, cv; cellset=getcellset(grid, "A"))
    ad2_B = AssemblyDomain("B", dh, m2, cv; cellset=getcellset(grid, "B"))
    b1, s1 = setup_assembly([ad1_A, ad1_B]; threading=true)
    b2, s2 = setup_assembly([ad2_A, ad2_B]; threading=true)
    K_test3 = test_modified_material(dh, b1, b2, s1, k1, k2)

    @test K_test1 ≈ K_test2
    @test K_test1 ≈ K_test3    
end