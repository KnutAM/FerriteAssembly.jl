@testset "Assemblers" begin
    # The Ferrite assembler has been tested other places, so we use that as the correct one here 
    grid = generate_grid(Quadrilateral, (10,8))
    ip = Lagrange{RefQuadrilateral,1}()
    dh = DofHandler(grid); add!(dh, :p, ip); add!(dh, :u, ip^2); close!(dh)
    ch = ConstraintHandler(dh); 
    add!(ch, Dirichlet(:p, getfacetset(grid, "left"), Returns(1.0))); 
    add!(ch, Dirichlet(:u, getfacetset(grid, "top"), Returns(Vec((2.0, 3.0)))))
    close!(ch)
    qr = QuadratureRule{RefQuadrilateral}(2); ip = Lagrange{RefQuadrilateral,1}()
    cv = (p=CellValues(qr, ip), u=CellValues(qr, ip^2))
    material = EE.PoroElasticPlaneStrain(;E=rand(), ν=rand()/2, k=rand(), α=rand(), β=rand())
    buffer = setup_domainbuffer(DomainSpec(dh, material, cv))

    set_time_increment!(buffer, 1.0)
    a = rand(ndofs(dh))
    aold = zeros(ndofs(dh))

    K0 = allocate_matrix(dh); r0 = zeros(ndofs(dh))
    ferrite_assembler = start_assemble(K0, r0)
    work!(ferrite_assembler, buffer; a=a, aold=aold)
    
    r1 = similar(r0)
    r_assembler = ReAssembler(r1)
    work!(r_assembler, buffer; a=a, aold=aold)

    K2 = similar(K0); r2 = similar(r0)
    @test_throws ArgumentError KeReAssembler(K2, r2; ch=ch) # apply_zero must be given
    kr_assembler = KeReAssembler(K2, r2)
    work!(kr_assembler, buffer; a=a, aold=aold)

    K3 = similar(K0); r3 = similar(r0)
    kr_assembler_ch0 = KeReAssembler(K3, r3; ch=ch, apply_zero=true)
    work!(kr_assembler_ch0, buffer; a=a, aold=aold)

    K4 = similar(K0); r4 = similar(r0)
    kr_assembler_ch1 = KeReAssembler(K4, r4; ch=ch, apply_zero=false)
    work!(kr_assembler_ch1, buffer; a=a, aold=aold)

    @test r1 ≈ r0
    @test K2 ≈ K0
    @test r2 ≈ r0
    
    # Note: Modification of stiffness different for local application, only solution can be checked. 
    K03 = copy(K0); r03 = copy(r0)
    apply_zero!(K03, r03, ch) 
    a03 = K03\r03
    a3 = K3\r3
    @test a03 ≈ a3

    K04 = copy(K0); r04 = copy(r0)
    apply!(K04, r04, ch)
    a04 = K04 \ r04
    a4 = K4\r4
    @test a04 ≈ a4

    @testset "can_thread with affine constraints" begin
        # Affine constraints make local `apply_assemble!` write to master dofs
        # outside the current cell, which is unsafe to thread over mesh coloring.
        struct OnesMaterial end
        function FerriteAssembly.element_routine!(Ke, re, state, ae, ::OnesMaterial, cv, buffer)
            fill!(Ke, 1.0)
            fill!(re, 1.0)
        end

        grid_line = generate_grid(Line, (8,))
        ip_line = Lagrange{RefLine, 1}()
        dh_line = close!(add!(DofHandler(grid_line), :u, ip_line))
        cv_line = CellValues(QuadratureRule{RefLine}(1), ip_line)

        ch_affine = ConstraintHandler(dh_line)
        for dof in 2:ndofs(dh_line)
            add!(ch_affine, AffineConstraint(dof, [1 => 1.0], 0.0))
        end
        close!(ch_affine)

        ch_dbc = ConstraintHandler(dh_line)
        add!(ch_dbc, Dirichlet(:u, Set([1]), Returns(0.0)))
        close!(ch_dbc)

        Kd = allocate_matrix(dh_line); rd = zeros(ndofs(dh_line))
        @test FA.can_thread(KeReAssembler(Kd, rd))
        @test FA.can_thread(KeReAssembler(Kd, rd; ch=ch_dbc, apply_zero=true))
        @test !FA.can_thread(KeReAssembler(Kd, rd; ch=ch_affine, apply_zero=true))

        # can_thread is cached at construction, so an open (not-yet-closed) `ch` must be
        # rejected: otherwise constraints added after construction wouldn't be reflected.
        ch_open = ConstraintHandler(dh_line)
        add!(ch_open, AffineConstraint(2, [1 => 1.0], 0.0))
        @test !Ferrite.isclosed(ch_open)
        @test_throws ArgumentError KeReAssembler(Kd, rd; ch=ch_open, apply_zero=true)

        db_sequential = setup_domainbuffer(DomainSpec(dh_line, OnesMaterial(), cv_line); threading=false)
        db_threaded = setup_domainbuffer(DomainSpec(dh_line, OnesMaterial(), cv_line); threading=true, num_tasks=4)

        Ka = allocate_matrix(dh_line, ch_affine); ra = zeros(ndofs(dh_line))
        asm_sequential = KeReAssembler(Ka, ra; ch=ch_affine, apply_zero=true)
        work!(asm_sequential, db_sequential)

        for _ in 1:5 # repeat to make a would-be race visible
            Kb = allocate_matrix(dh_line, ch_affine); rb = zeros(ndofs(dh_line))
            asm_threaded = KeReAssembler(Kb, rb; ch=ch_affine, apply_zero=true)
            work!(asm_threaded, db_threaded)
            @test Kb ≈ Ka
            @test rb ≈ ra
        end
    end
end