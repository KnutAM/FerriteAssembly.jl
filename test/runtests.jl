using Ferrite, FerriteAssembly
using LinearAlgebra, SparseArrays
using Test
import FerriteAssembly as FA
import FerriteAssembly.ExampleElements as EE
import MaterialModelsBase as MMB

include("states.jl") 
include("threading_utils.jl")
include("heatequation.jl")
include("example_elements.jl")
include("integration.jl")
include("scaling.jl")
include("CustomStiffness.jl")


@testset "Errors" begin
    printstyled("=== Testing will give expected error messages, ok if tests pass! ===\n"; color=:green, bold=true)
    grid = generate_grid(Quadrilateral, (2,2))
    dh = DofHandler(grid); add!(dh, :u, 1); close!(dh)
    cv = CellScalarValues(QuadratureRule{2, RefCube}(2), Lagrange{2, RefCube, 1}());
    K = create_sparsity_pattern(dh)
    r = zeros(ndofs(dh))
    a = zeros(ndofs(dh))
    # Check error when element_residual! is not defined. material=nothing
    buffer, new_states, old_states = setup_assembly(dh, nothing, cv)
    buffer_ad, new_states, old_states = setup_assembly(dh, nothing, cv; autodiffbuffer=true)
    
    exception = ErrorException("Did not find correctly defined element_routine! or element_residual!")
    assembler = start_assemble(K, r)
    @test_throws exception doassemble!(assembler, new_states, buffer; a=a)
    @test_throws exception doassemble!(assembler, new_states, buffer_ad; a=a)
    
    # Check error when trying to insert a dual number into the state variables
    # Note: Not optimal as we should actually check the printed message
    struct _TestMaterial end 
    FerriteAssembly.create_cell_state(::_TestMaterial, args...) = zeros(1)
    
    function FerriteAssembly.element_residual!(re, state, ae, m::_TestMaterial, args...)
        state[1] = first(ae) # Not allowed, must be state[1] = ForwardDiff.value(first(ae))
    end
    buffer_dualissue, old_states_dualissue, new_states_dualissue = setup_assembly(dh, _TestMaterial(), cv; autodiffbuffer=true)
    @test_throws MethodError doassemble!(assembler, new_states_dualissue, buffer_dualissue; a=a, old_states=old_states_dualissue)

    printstyled("================== End of expected error messages ==================\n"; color=:green, bold=true)
end

@testset "Miscellaneous" begin
    grid = generate_grid(Quadrilateral, (2,2))
    dh = DofHandler(grid); add!(dh, :u, 1); close!(dh)
    cv = CellScalarValues(QuadratureRule{2, RefCube}(2), Lagrange{2, RefCube, 1}());
    @testset "get_material" begin
        material = zeros(1) #dummy
        buffer, new_states, old_states = setup_assembly(dh, material, cv; threading=Val(false))
        buffer_threaded, _, _ = setup_assembly(dh, material, cv; threading=Val(true))
        buffers, _ = setup_assembly([AssemblyDomain("a", dh, material, cv)]; threading=Val(false))
        @test material === FerriteAssembly.get_material(buffer)
        @test material === FerriteAssembly.get_material(buffer_threaded)
        @test material === FerriteAssembly.get_material(buffers, "a")
    end
end

# Print show warning at the end if running tests single-threaded. 
# During unit testing, exclude tests that must be multithreaded to pass.
Threads.nthreads() == 1 && @warn("Threads.nthreads() == 1: Run multithreaded for full test coverage")
# =#