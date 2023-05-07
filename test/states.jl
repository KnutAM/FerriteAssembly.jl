module TestStateModule
    using Ferrite, FerriteAssembly
    struct MatA end
    struct StateA # bitstype
        cellnr::Int
        quadnr::Int
    end
    FerriteAssembly.create_cell_state(::MatA, cv, args...) = [StateA(-1, 0) for _ in 1:getnquadpoints(cv)]
    function FerriteAssembly.element_residual!(re, states_new::Vector{StateA}, ae, ::MatA, cv, buffer)
        cellnr = cellid(buffer)
        for i in 1:getnquadpoints(cv)
            states_new[i] = StateA(cellnr, i)
        end
    end

    struct MatB{dim} end
    mutable struct StateB{dim} # not bitstype
        cellnr::Int
        const quad_coordinates::Vector{Vec{dim,Float64}}
    end
    FerriteAssembly.create_cell_state(::MatB{dim}, cv, args...) where dim = StateB(-1, [zero(Vec{dim}) for i in 1:getnquadpoints(cv)])
    function FerriteAssembly.element_residual!(re, states_new::StateB, ae, ::MatB, cv, buffer)
        coords = getcoordinates(buffer)
        states_new.cellnr = cellid(buffer)
        for i in 1:getnquadpoints(cv)
            x = spatial_coordinate(cv, i, coords)
            states_new.quad_coordinates[i] = x
        end
    end

    struct MatC end
    struct StateC
        counter::Int
    end
    FerriteAssembly.create_cell_state(::MatC, cv, args...) = [StateC(0) for _ in 1:getnquadpoints(cv)]
    function FerriteAssembly.element_residual!(re, states_new::Vector{StateC}, ae, ::MatC, cv, buffer)
        states_old = FerriteAssembly.get_state_old(buffer)
        for i in 1:getnquadpoints(cv)
            states_new[i] = StateC(states_old[i].counter + 1)
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
        a = copy(r)
        ip = Ferrite.default_interpolation(getcelltype(grid))
        RefShape = Ferrite.getrefshape(ip)
        cv = CellVectorValues(QuadratureRule{Dim,RefShape}(2), ip);

        # MatA: 
        # - Check correct values before and after update
        # - Check unaliased old and new after update_states!
        buffer, states_old, states_new = setup_assembly(dh, MatA(), cv)
        @test isa(states_old, Dict{Int,Vector{StateA}})
        @test states_old == states_new
        @test states_old[1] == [StateA(-1, 0) for _ in 1:getnquadpoints(cv)]
        doassemble!(r, states_new, states_old, buffer)
        @test states_old[1] == [StateA(-1, 0) for _ in 1:getnquadpoints(cv)] # Unchanged
        for cellnr in 1:getncells(grid)
            @test states_new[cellnr] == [StateA(cellnr, i) for i in 1:getnquadpoints(cv)]  # Updated
        end
        FerriteAssembly.update_states!(states_old, states_new)
        @test states_old == states_new          # Correctly updated values
        states_new[1][1] = StateA(0,0)
        @test states_old[1][1] == StateA(1,1)   # But not aliased
        allocs = @allocated FerriteAssembly.update_states!(states_old, states_new)
        @test allocs == 0 # Vector{T} where isbitstype(T) should not allocate (MatA fulfills this)

        # MatB (not bitstype)
        # - Check correct values before and after update
        # - Check unaliased old and new after update_states!
        buffer, states_old, states_new = setup_assembly(dh, MatB{Dim}(), cv)
        @test isa(states_old, Dict{Int,StateB{Dim}})
        @test states_old == states_new
        @test states_old[1] == StateB(-1, [zero(Vec{Dim}) for i in 1:getnquadpoints(cv)])
        doassemble!(K, r, states_new, states_old, buffer)
        @test states_old[1] == StateB(-1, [zero(Vec{Dim}) for i in 1:getnquadpoints(cv)]) # Unchanged
        for cellnr in 1:getncells(grid)
            coords = getcoordinates(grid, cellnr)
            x_values = [spatial_coordinate(cv, i, coords) for i in 1:getnquadpoints(cv)]
            @test states_new[cellnr] == StateB(cellnr, x_values)                          # Updated
        end
        FerriteAssembly.update_states!(states_old, states_new)
        @test states_old == states_new                  # Correctly updated values
        cellnr = rand(1:getncells(grid))
        coords = getcoordinates(grid, cellnr)
        x_values = [spatial_coordinate(cv, i, coords) for i in 1:getnquadpoints(cv)]
        states_new[cellnr] = StateB(0, -x_values)
        @test states_old[cellnr] == StateB(cellnr, x_values)   # But not aliased
        allocs = @allocated FerriteAssembly.update_states!(states_old, states_new)
        @test allocs > 0 # Vector{T} where !isbitstype(T) is expected to allocate. 
        # If this fails after improvements, that is just good, but then the docstring should be updated. 

        # MatC (accumulation), using threading as well
        colors = create_coloring(grid)
        buffer, states_old, states_new = setup_assembly(dh, MatC(), cv; colors=colors)
        @test isa(states_old, Dict{Int,Vector{StateC}})
        @test states_old == states_new
        @test states_old[1][1] == StateC(0)
        doassemble!(K, r, states_new, states_old, buffer)
        @test states_old[1][1] == StateC(0)
        @test states_new[1][1] == StateC(1) # Added 1
        FerriteAssembly.update_states!(states_old, states_new)
        @test states_old[1][1] == StateC(1) # Updated
        states_new[1][1] = StateC(0)        # Set states_new to zero to test aliasing and update in next assembly
        @test states_old[1][1] == StateC(1) # But not aliased
        doassemble!(K, r, states_new, states_old, buffer)
        @test states_new[1][1] == StateC(2) # Added 1 from states_old[1][1] (not from states_new[1][1] which was StateC(0))
        for cellnr in 1:getncells(grid)
            @test states_new[cellnr][2] == StateC(2) # Check that all are updated
        end
        allocs = @allocated FerriteAssembly.update_states!(states_old, states_new)
        @test allocs == 0 # Vector{T} where isbitstype(T) should not allocate (MatC fulfills this)
    end
end
