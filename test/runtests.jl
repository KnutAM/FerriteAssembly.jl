using Ferrite, FerriteAssembly
using SparseArrays
using Test

@testset "FerriteAssembly.jl" begin
    include("heatequation.jl")

    @testset "Errors" begin
        grid = generate_grid(Quadrilateral, (2,2))
        dh = DofHandler(grid); push!(dh, :u, 1); close!(dh)
        cv = CellScalarValues(QuadratureRule{2, RefCube}(2), Lagrange{2, RefCube, 1}());
        states = create_states(dh)
        K = create_sparsity_pattern(dh)
        r = zeros(ndofs(dh))
        a = zeros(ndofs(dh))
        # Check error when element_residual! is not defined. 
        # Note: Not optimal as we should actually check the printed message
        cb = CellBuffer(dh, cv, nothing)    # material=nothing not supported
        cb_ad = AutoDiffCellBuffer(states, dh, cv, nothing)
        exception = ErrorException("Did not find correctly defined element_routine! or element_residual!")
        @test_throws exception doassemble!(start_assemble(K,r), cb, states, dh, a)
        @test_throws exception doassemble!(start_assemble(K,r), cb_ad, states, dh, a)
        
        # Check error when trying to insert a dual number into the state variables
        # Note: Not optimal as we should actually check the printed message
        states_dualissue = create_states(dh, Returns(0.0), cv)
        struct _TestMaterial end 
        
        function FerriteAssembly.element_residual!(re, state, ae, m::_TestMaterial, args...)
            state[1] = first(ae) # Not allowed, must be state[1] = ForwardDiff.value(first(ae))
        end
        cb_dualissue = AutoDiffCellBuffer(states_dualissue, dh, cv, _TestMaterial())
        @test_throws MethodError doassemble!(start_assemble(K,r), cb_dualissue, states_dualissue, dh, a) 

        # Check that error thrown for mismatch between parallel buffer sizes and number of threads 
        @test Threads.nthreads() > 1    # OK to fail locally, but for tests to be proper, we need to run CI multithreaded
        if Threads.nthreads() > 1 # The following test will fail unless running more than one thread
            cellbuffers = create_threaded_CellBuffers(cb)
            assemblers = create_threaded_assemblers(K, r)
            cellbuffers_n1 = create_threaded_CellBuffers(cb; nthreads=1)
            assemblers_n1 = create_threaded_assemblers(K, r; nthreads=1)
            colors = create_coloring(dh.grid)
            for (as,cb) in ((assemblers_n1, cellbuffers_n1), (assemblers_n1, cellbuffers), (assemblers, cellbuffers_n1))
                @test_throws DimensionMismatch doassemble!(as, cb, states, dh, colors, a)
            end
        end
    end

end
