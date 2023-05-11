module TestStateModule
    using Ferrite, FerriteAssembly
    struct MatA end
    struct StateA # bitstype
        cellnr::Int
        quadnr::Int
    end
    FerriteAssembly.create_cell_state(::MatA, cv, args...) = [StateA(-1, 0) for _ in 1:getnquadpoints(cv)]
    function FerriteAssembly.element_residual!(re, new_states::Vector{StateA}, ae, ::MatA, cv, buffer)
        cellnr = cellid(buffer)
        for i in 1:getnquadpoints(cv)
            new_states[i] = StateA(cellnr, i)
        end
    end

    struct MatB{dim} end
    mutable struct StateB{dim} # not bitstype
        cellnr::Int
        const quad_coordinates::Vector{Vec{dim,Float64}}
    end
    FerriteAssembly.create_cell_state(::MatB{dim}, cv, args...) where dim = StateB(-1, [zero(Vec{dim}) for i in 1:getnquadpoints(cv)])
    function FerriteAssembly.element_residual!(re, new_states::StateB, ae, ::MatB, cv, buffer)
        coords = getcoordinates(buffer)
        new_states.cellnr = cellid(buffer)
        for i in 1:getnquadpoints(cv)
            x = spatial_coordinate(cv, i, coords)
            new_states.quad_coordinates[i] = x
        end
    end

    struct MatC end
    struct StateC
        counter::Int
    end
    FerriteAssembly.create_cell_state(::MatC, cv, args...) = [StateC(0) for _ in 1:getnquadpoints(cv)]
    function FerriteAssembly.element_residual!(re, new_states::Vector{StateC}, ae, ::MatC, cv, buffer)
        old_states = FerriteAssembly.get_old_state(buffer)
        for i in 1:getnquadpoints(cv)
            new_states[i] = StateC(old_states[i].counter + 1)
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
        dh = DofHandler(grid); add!(dh, :u, Dim); close!(dh);

        K = create_sparsity_pattern(dh)
        r = zeros(ndofs(dh));
        kr_assembler = start_assemble(K, r)
        r_assembler = FerriteAssembly.ReAssembler(r)
        a = copy(r)
        ip = Ferrite.default_interpolation(getcelltype(grid))
        RefShape = Ferrite.getrefshape(ip)
        cv = CellVectorValues(QuadratureRule{Dim,RefShape}(2), ip);

        # MatA: 
        # - Check correct values before and after update
        # - Check unaliased old and new after update_states!
        buffer, old_states, new_states = setup_assembly(dh, MatA(), cv)
        @test isa(old_states, Dict{Int,Vector{StateA}})
        @test old_states == new_states
        @test old_states[1] == [StateA(-1, 0) for _ in 1:getnquadpoints(cv)]
        doassemble!(r_assembler, new_states, buffer; old_states=old_states)
        @test old_states[1] == [StateA(-1, 0) for _ in 1:getnquadpoints(cv)] # Unchanged
        for cellnr in 1:getncells(grid)
            @test new_states[cellnr] == [StateA(cellnr, i) for i in 1:getnquadpoints(cv)]  # Updated
        end
        update_states!(old_states, new_states)
        @test old_states == new_states          # Correctly updated values
        new_states[1][1] = StateA(0,0)
        @test old_states[1][1] == StateA(1,1)   # But not aliased
        allocs = @allocated update_states!(old_states, new_states)
        @test allocs == 0 # Vector{T} where isbitstype(T) should not allocate (MatA fulfills this)

        # MatB (not bitstype)
        # - Check correct values before and after update
        # - Check unaliased old and new after update_states!
        buffer, old_states, new_states = setup_assembly(dh, MatB{Dim}(), cv)
        @test isa(old_states, Dict{Int,StateB{Dim}})
        @test old_states == new_states
        @test old_states[1] == StateB(-1, [zero(Vec{Dim}) for i in 1:getnquadpoints(cv)])
        doassemble!(kr_assembler, new_states, buffer; old_states=old_states)
        @test old_states[1] == StateB(-1, [zero(Vec{Dim}) for i in 1:getnquadpoints(cv)]) # Unchanged
        for cellnr in 1:getncells(grid)
            coords = getcoordinates(grid, cellnr)
            x_values = [spatial_coordinate(cv, i, coords) for i in 1:getnquadpoints(cv)]
            @test new_states[cellnr] == StateB(cellnr, x_values)                          # Updated
        end
        update_states!(old_states, new_states)
        @test old_states == new_states                  # Correctly updated values
        cellnr = rand(1:getncells(grid))
        coords = getcoordinates(grid, cellnr)
        x_values = [spatial_coordinate(cv, i, coords) for i in 1:getnquadpoints(cv)]
        new_states[cellnr] = StateB(0, -x_values)
        @test old_states[cellnr] == StateB(cellnr, x_values)   # But not aliased
        allocs = @allocated update_states!(old_states, new_states)
        @test allocs > 0 # Vector{T} where !isbitstype(T) is expected to allocate. 
        # If this fails after improvements, that is just good, but then the docstring should be updated. 

        # MatC (accumulation), using threading as well
        colors = create_coloring(grid)
        buffer, old_states, new_states = setup_assembly(dh, MatC(), cv; colors=colors)
        @test isa(old_states, Dict{Int,Vector{StateC}})
        @test old_states == new_states
        @test old_states[1][1] == StateC(0)
        doassemble!(kr_assembler, new_states, buffer; old_states=old_states)
        @test old_states[1][1] == StateC(0)
        @test new_states[1][1] == StateC(1) # Added 1
        update_states!(old_states, new_states)
        @test old_states[1][1] == StateC(1) # Updated
        new_states[1][1] = StateC(0)        # Set new_states to zero to test aliasing and update in next assembly
        @test old_states[1][1] == StateC(1) # But not aliased
        doassemble!(kr_assembler, new_states, buffer; old_states=old_states)
        @test new_states[1][1] == StateC(2) # Added 1 from old_states[1][1] (not from new_states[1][1] which was StateC(0))
        for cellnr in 1:getncells(grid)
            @test new_states[cellnr][2] == StateC(2) # Check that all are updated
        end
        allocs = @allocated update_states!(old_states, new_states)
        @test allocs == 0 # Vector{T} where isbitstype(T) should not allocate (MatC fulfills this)
    end

    # Smoke-test of update_states! for nothing states (and check no allocations)
    snew = Dict(i=>nothing for i in 1:10)
    sold = deepcopy(snew)
    update_states!(sold, snew) # Compile
    allocs = @allocated update_states!(sold, snew)
    @test allocs == 0
    snew = Dict(string(k)=>Dict(i=>nothing for i in k:2:10) for k in 1:2)
    sold = deepcopy(snew)
    update_states!(sold, snew) # Compile
    allocs = @allocated update_states!(sold, snew)
    @test allocs == 0
end
