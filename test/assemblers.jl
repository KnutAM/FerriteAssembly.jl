@testset "Assemblers" begin
    # The Ferrite assembler has been tested other places, so we use that as the correct one here 
    grid = generate_grid(Quadrilateral, (10,8))
    ip = Lagrange{RefQuadrilateral,1}()
    dh = DofHandler(grid); add!(dh, :p, ip); add!(dh, :u, ip^2); close!(dh)
    ch = ConstraintHandler(dh); 
    add!(ch, Dirichlet(:p, getfaceset(grid, "left"), Returns(1.0))); 
    add!(ch, Dirichlet(:u, getfaceset(grid, "top"), Returns(Vec((2.0, 3.0)))))
    close!(ch)
    qr = QuadratureRule{RefQuadrilateral}(2); ip = Lagrange{RefQuadrilateral,1}()
    cv = (p=CellValues(qr, ip), u=CellValues(qr, ip^2))
    material = EE.PoroElasticPlaneStrain(;E=rand(), ν=rand()/2, k=rand(), α=rand(), β=rand())
    buffer = setup_domainbuffer(DomainSpec(dh, material, cv))

    set_time_increment!(buffer, 1.0)
    a = rand(ndofs(dh))
    aold = zeros(ndofs(dh))

    K0 = create_sparsity_pattern(dh); r0 = zeros(ndofs(dh))
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
end