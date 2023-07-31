@testset "Errors" begin
    printstyled("=== Testing will give expected error messages, ok if tests pass! ===\n"; color=:green, bold=true)
    grid = generate_grid(Quadrilateral, (2,2))
    dh = DofHandler(grid); add!(dh, :u, 1); close!(dh)
    cv = CellScalarValues(QuadratureRule{2, RefCube}(2), Lagrange{2, RefCube, 1}());
    K = create_sparsity_pattern(dh)
    r = zeros(ndofs(dh))
    a = zeros(ndofs(dh))
    # Check error when element_residual! is not defined. material=nothing
    grid_domain = DomainSpec(dh, nothing, cv)
    buffer = setup_domainbuffer(grid_domain)
    buffer_ad = setup_domainbuffer(grid_domain; autodiffbuffer=true)
    
    exception = ErrorException("Did not find correctly defined element_routine! or element_residual!")
    assembler = start_assemble(K, r)
    @test_throws exception work!(assembler, buffer; a=a)
    @test_throws exception work!(assembler, buffer_ad; a=a)
    
    # Check error when trying to insert a dual number into the state variables
    # Note: Not optimal as we should actually check the printed message
    struct _TestMaterial end 
    FerriteAssembly.create_cell_state(::_TestMaterial, args...) = zeros(1)
    
    function FerriteAssembly.element_residual!(re, state, ae, m::_TestMaterial, args...)
        state[1] = first(ae) # Not allowed, must be state[1] = ForwardDiff.value(first(ae))
    end
    buffer_dualissue = setup_domainbuffer(DomainSpec(dh, _TestMaterial(), cv); autodiffbuffer=true)
    @test_throws MethodError work!(assembler, buffer_dualissue; a=a)

    printstyled("================== End of expected error messages ==================\n"; color=:green, bold=true)
end
