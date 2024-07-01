module TestStateModule
    using Ferrite, FerriteAssembly
    struct MatA end
    struct StateA # bitstype
        cellnr::Int
        quadnr::Int
    end
    FerriteAssembly.create_cell_state(::MatA, cv, args...) = [StateA(-1, 0) for _ in 1:getnquadpoints(cv)]
    function FerriteAssembly.element_residual!(re, states::Vector{StateA}, ae, ::MatA, cv, buffer)
        cellnr = cellid(buffer)
        for i in 1:getnquadpoints(cv)
            states[i] = StateA(cellnr, i)
        end
    end

    struct MatB{dim} end
    mutable struct StateB{dim} # not bitstype
        cellnr::Int
        const quad_coordinates::Vector{Vec{dim,Float64}}
    end
    FerriteAssembly.create_cell_state(::MatB{dim}, cv, args...) where dim = StateB(-1, [zero(Vec{dim}) for i in 1:getnquadpoints(cv)])
    function FerriteAssembly.element_residual!(re, states::StateB, ae, ::MatB, cv, buffer)
        coords = getcoordinates(buffer)
        states.cellnr = cellid(buffer)
        for i in 1:getnquadpoints(cv)
            x = spatial_coordinate(cv, i, coords)
            states.quad_coordinates[i] = x
        end
    end

    struct MatC end
    struct StateC
        counter::Int
    end
    FerriteAssembly.create_cell_state(::MatC, cv, args...) = [StateC(0) for _ in 1:getnquadpoints(cv)]
    function FerriteAssembly.element_residual!(re, states::Vector{StateC}, ae, ::MatC, cv, buffer)
        old_states = FerriteAssembly.get_old_state(buffer)
        for i in 1:getnquadpoints(cv)
            states[i] = StateC(old_states[i].counter + 1)
        end
    end

    # With contained vectors, comparison gives false in all cases even with equal values...
    Base.:(==)(a::StateB, b::StateB) = (a.cellnr==b.cellnr && (mapreduce((ax, bx)->ax==bx, *, a.quad_coordinates, b.quad_coordinates)))

end

@testset "state variables" begin
    # Defs
    import .TestStateModule: MatA, MatB, MatC, StateA, StateB, StateC
    
    for (CT, Dim) in ((Line, 1), (QuadraticTriangle, 2), (Hexahedron, 3))
        grid = generate_grid(CT, ntuple(_->3, Dim))
        ip = geometric_interpolation(CT)
        dh = DofHandler(grid); add!(dh, :u, ip^Dim); close!(dh);

        K = allocate_matrix(dh)
        r = zeros(ndofs(dh));
        kr_assembler = start_assemble(K, r)
        r_assembler = FerriteAssembly.ReAssembler(r)
        a = copy(r)
        RefShape = Ferrite.getrefshape(ip)
        cv = CellValues(QuadratureRule{RefShape}(2), ip^Dim, ip);

        # MatA: 
        # - Check correct values before and after update
        # - Check unaliased old and new after update_states!
        buffer = setup_domainbuffer(DomainSpec(dh, MatA(), cv))
        states = FerriteAssembly.get_state(buffer)
        old_states = FerriteAssembly.get_old_state(buffer)
        @test isa(old_states, Dict{Int,Vector{StateA}})
        @test old_states == states
        @test old_states[1] == [StateA(-1, 0) for _ in 1:getnquadpoints(cv)]
        work!(r_assembler, buffer)
        @test old_states[1] == [StateA(-1, 0) for _ in 1:getnquadpoints(cv)] # Unchanged
        for cellnr in 1:getncells(grid)
            @test states[cellnr] == [StateA(cellnr, i) for i in 1:getnquadpoints(cv)]  # Updated
        end
        update_states!(buffer)
        @test old_states == states          # Correctly updated values
        states[1][1] = StateA(0,0)
        @test old_states[1][1] == StateA(1,1)   # But not aliased
        allocs = @allocated update_states!(buffer)
        @test allocs == 0 # Vector{T} where isbitstype(T) should not allocate (MatA fulfills this)

        # MatB (not bitstype)
        # - Check correct values before and after update
        # - Check unaliased old and new after update_states!
        buffer = setup_domainbuffer(DomainSpec(dh, MatB{Dim}(), cv))
        states = FerriteAssembly.get_state(buffer)
        old_states = FerriteAssembly.get_old_state(buffer)
        @test isa(old_states, Dict{Int,StateB{Dim}})
        @test old_states == states
        @test old_states[1] == StateB(-1, [zero(Vec{Dim}) for i in 1:getnquadpoints(cv)])
        work!(kr_assembler, buffer)
        @test old_states[1] == StateB(-1, [zero(Vec{Dim}) for i in 1:getnquadpoints(cv)]) # Unchanged
        for cellnr in 1:getncells(grid)
            coords = getcoordinates(grid, cellnr)
            x_values = [spatial_coordinate(cv, i, coords) for i in 1:getnquadpoints(cv)]
            @test states[cellnr] == StateB(cellnr, x_values)                          # Updated
        end
        update_states!(buffer)
        @test old_states == states                  # Correctly updated values
        cellnr = rand(1:getncells(grid))
        coords = getcoordinates(grid, cellnr)
        x_values = [spatial_coordinate(cv, i, coords) for i in 1:getnquadpoints(cv)]
        states[cellnr] = StateB(0, -x_values)
        @test old_states[cellnr] == StateB(cellnr, x_values)   # But not aliased
        allocs = @allocated update_states!(buffer)
        @test allocs > 0 # Vector{T} where !isbitstype(T) is expected to allocate. 
        # If this fails after improvements, that is just good, but then the docstring should be updated. 

        # MatC (accumulation), using threading as well
        colors = create_coloring(grid)
        buffer = setup_domainbuffer(DomainSpec(dh, MatC(), cv; colors=colors))
        states = FerriteAssembly.get_state(buffer)
        old_states = FerriteAssembly.get_old_state(buffer)
        @test isa(old_states, Dict{Int,Vector{StateC}})
        @test old_states == states
        @test old_states[1][1] == StateC(0)
        work!(kr_assembler, buffer)
        @test old_states[1][1] == StateC(0)
        @test states[1][1] == StateC(1) # Added 1
        update_states!(buffer)
        @test old_states[1][1] == StateC(1) # Updated
        states[1][1] = StateC(0)        # Set states to zero to test aliasing and update in next assembly
        @test old_states[1][1] == StateC(1) # But not aliased
        work!(kr_assembler, buffer)
        @test states[1][1] == StateC(2) # Added 1 from old_states[1][1] (not from states[1][1] which was StateC(0))
        for cellnr in 1:getncells(grid)
            @test states[cellnr][2] == StateC(2) # Check that all are updated
        end
        allocs = @allocated update_states!(buffer)
        @test allocs == 0 # Vector{T} where isbitstype(T) should not allocate (MatC fulfills this)
    end

    ip = Lagrange{RefTriangle,1}()
    dh = DofHandler(generate_grid(Triangle, (2,2))); add!(dh, :u, ip); close!(dh);
    
    # Smoke-test of update_states! for nothing states (and check no allocations)
    cv = CellValues(QuadratureRule{RefTriangle}(2), ip)
    buffer = setup_domainbuffer(DomainSpec(dh, nothing, cv))
    @test isa(FerriteAssembly.get_state(buffer), Dict{Int,Nothing})
    update_states!(buffer) # Compile
    allocs = @allocated update_states!(buffer)
    @test allocs == 0

    gda = DomainSpec(dh, nothing, cv; set=1:getncells(dh.grid)÷2)
    gdb = DomainSpec(dh, nothing, cv; set=setdiff!(Set(1:getncells(dh.grid)), gda.set))
    buffers = setup_domainbuffers(Dict("a"=>gda, "b"=>gdb))
    update_states!(buffers) # Compile
    allocs = @allocated update_states!(buffers)
    @test allocs == 0
end
