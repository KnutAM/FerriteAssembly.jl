using Ferrite, FerriteAssembly
using ForwardDiff
using SparseArrays
using Test
import FerriteAssembly as FA
import FerriteAssembly.ExampleElements as EE
import MaterialModelsBase as MMB
using Logging

include("replacements.jl")
include("states.jl") 
include("threading_utils.jl")
include("heatequation.jl")
include("example_elements.jl")
include("setup.jl")
include("assemblers.jl")
include("facet_assembly.jl")
include("integration.jl")
include("quadpoint_evaluation.jl")
include("load_handler.jl")
include("scaling.jl")
include("errors.jl")

@testset "Miscellaneous" begin
    grid = generate_grid(Quadrilateral, (2,2))
    ip = Lagrange{RefQuadrilateral, 1}()
    dh = DofHandler(grid); add!(dh, :u, ip); close!(dh)
    cv = CellValues(QuadratureRule{RefQuadrilateral}(2), ip);
    @testset "get functions" begin
        material = zeros(1) #dummy
        grid_domain = DomainSpec(dh, material, cv)
        buffer = setup_domainbuffer(grid_domain; threading=false)
        buffer_threaded = setup_domainbuffer(grid_domain; threading=true)
        buffers = setup_domainbuffers(Dict("a"=>grid_domain); threading=false)
        @test material === FerriteAssembly.get_material(buffer)
        @test material === FerriteAssembly.get_material(buffer_threaded)
        @test material === FerriteAssembly.get_material(buffers, "a")

        @test FerriteAssembly.get_itembuffer(buffer) isa FerriteAssembly.CellBuffer
        @test FerriteAssembly.get_itembuffer(buffer_threaded) isa FerriteAssembly.TaskLocals{<:FerriteAssembly.CellBuffer}
        @test FerriteAssembly.get_itembuffer(buffers, "a") isa FerriteAssembly.CellBuffer

        @test FerriteAssembly.get_dofhandler(buffer) === dh
        @test FerriteAssembly.get_dofhandler(buffer_threaded) === dh
        @test FerriteAssembly.get_dofhandler(buffers) === dh

        @test FerriteAssembly.get_state(buffer, 1) === nothing
        @test FerriteAssembly.get_old_state(buffer, 1) === nothing
        
    end

    @testset "utility functions" begin
        x = rand()
        xd = ForwardDiff.Dual(x, rand(3)...)
        @test x === FerriteAssembly.remove_dual(x)
        @test x === FerriteAssembly.remove_dual(xd)

        t = rand(Tensor{2,3})
        td = Tensor{2,3}((i, j) -> ForwardDiff.Dual(t[i, j], rand(), rand()))
        @test t === FerriteAssembly.remove_dual(t)
        @test t === FerriteAssembly.remove_dual(td)
    end
end

# Print show warning at the end if running tests single-threaded. 
# During unit testing, exclude tests that must be multithreaded to pass.
Threads.nthreads() == 1 && @warn("Threads.nthreads() == 1: Run multithreaded for full test coverage")
# =#