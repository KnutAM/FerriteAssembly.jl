using Ferrite, FerriteAssembly
using SparseArrays
using Test
import FerriteAssembly as FA
import FerriteAssembly.ExampleElements as EE
import MaterialModelsBase as MMB

include("states.jl") 
include("heatequation.jl")
include("example_elements.jl")


@testset "Errors" begin
    printstyled("=== Testing will give expected error messages, ok if tests pass! ===\n"; color=:green, bold=true)
    grid = generate_grid(Quadrilateral, (2,2))
    dh = DofHandler(grid); add!(dh, :u, 1); close!(dh)
    cv = CellScalarValues(QuadratureRule{2, RefCube}(2), Lagrange{2, RefCube, 1}());
    K = create_sparsity_pattern(dh)
    r = zeros(ndofs(dh))
    a = zeros(ndofs(dh))
    # Check error when element_residual! is not defined. material=nothing
    buffer, states_old, states_new = setup_assembly(dh, nothing, cv)
    buffer_ad, states_old, states_new = setup_assembly(dh, nothing, cv; autodiffbuffer=true)
    
    exception = ErrorException("Did not find correctly defined element_routine! or element_residual!")
    @test_throws exception doassemble!(K, r, states_new, states_old, buffer; a=a)
    @test_throws exception doassemble!(K, r, states_new, states_old, buffer_ad; a=a)
    
    # Check error when trying to insert a dual number into the state variables
    # Note: Not optimal as we should actually check the printed message
    struct _TestMaterial end 
    FerriteAssembly.create_cell_state(::_TestMaterial, args...) = zeros(1)
    
    function FerriteAssembly.element_residual!(re, state, ae, m::_TestMaterial, args...)
        state[1] = first(ae) # Not allowed, must be state[1] = ForwardDiff.value(first(ae))
    end
    buffer_dualissue, states_old_dualissue, states_new_dualissue = setup_assembly(dh, _TestMaterial(), cv; autodiffbuffer=true)
    @test_throws MethodError doassemble!(K, r, states_old_dualissue, states_new_dualissue, buffer_dualissue; a=a) 

    printstyled("================== End of expected error messages ==================\n"; color=:green, bold=true)
end

# Print show warning at the end if running tests single-threaded. 
# During unit testing, exclude tests that must be multithreaded to pass.
Threads.nthreads() == 1 && @warn("Threads.nthreads() == 1: Run multithreaded for full test coverage")