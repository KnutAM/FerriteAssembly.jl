@testset "set_chunks" begin
    for approx_num_points in (10,100,1000)
        for num_tasks in (1,2,4,8,16,32)
            set = unique!([rand(1:(approx_num_points*10)) for _ in 1:approx_num_points])
            merged_chunks = Set{Int}()
            for set_chunk in FerriteAssembly.split_in_chunks(set; num_tasks=num_tasks)
                union!(merged_chunks, set_chunk)
            end
            @test merged_chunks == Set(set)
        end
    end
end

@testset "get_chunk with empty chunks, PR89" begin
    # A empty chunk (e.g. supplied by a user via `chunks=...`, or
    # produced by an unlucky coloring/splitting) could previously be mistaken by a
    # worker as queue exhaustion. Check that this doesn't occur.
    for num_tasks in (1, 2, 4)
        # `num_tasks` consecutive leading empty chunks guarantee that, under the
        # old (buggy) implementation where an empty chunk was indistinguishable
        # from exhaustion, every one of the `num_tasks` workers would see an
        # empty chunk and exit before ever reaching the real work below.
        leading_empty = [Int[] for _ in 1:num_tasks]
        chunk_vector = [leading_empty..., [1, 2], Int[], [3], Int[]]
        taskchunks = FerriteAssembly.TaskChunks(chunk_vector)
        visited = Int[]
        # Bound the number of retrievals instead of looping until `nothing`:
        # if a regression reintroduces `Int[]` as the exhaustion sentinel,
        # this must fail promptly rather than hang forever.
        for _ in 1:length(chunk_vector)
            taskchunk = FerriteAssembly.get_chunk(taskchunks)
            @test taskchunk !== nothing
            taskchunk === nothing && break
            append!(visited, taskchunk)
        end
        @test sort(visited) == [1, 2, 3]
        # Queue must remain exhausted (not silently reset) once drained.
        @test FerriteAssembly.get_chunk(taskchunks) === nothing
    end
end

@testset "nonpositive num_tasks rejected (BUG-007)" begin
    # `num_tasks <= 0` previously created empty task-local arrays, silently
    # skipping all work instead of raising an error at setup.
    grid = generate_grid(Quadrilateral, (2, 1))
    ip = Lagrange{RefQuadrilateral, 1}()
    dh = close!(add!(DofHandler(grid), :u, ip))
    cv = CellValues(QuadratureRule{RefQuadrilateral}(2), ip)
    for num_tasks in (0, -1)
        @test_throws ArgumentError setup_domainbuffer(
            DomainSpec(dh, nothing, cv); threading=true, num_tasks=num_tasks)
    end
    for num_tasks in (1, Threads.nthreads() + 10)
        db = setup_domainbuffer(
            DomainSpec(dh, nothing, cv); threading=true, num_tasks=num_tasks)
        ig = SimpleIntegrator(Returns(1.0), 0.0)
        work!(ig, db)
        @test ig.val ≈ 4.0 # total area of the two-cell grid
    end
end

# `create_chunks` previously dispatched on concrete `Ferrite.Grid`, so any other
# `Ferrite.AbstractGrid` implementation (e.g. FerriteIGA.jl's `BezierGrid`) could not
# be used with `threading=true`, even though `create_coloring` and the generic
# `AbstractGrid` accessors it relies on only require `.cells` and `.nodes` fields.
struct DummyGrid{C,N} <: Ferrite.AbstractGrid{2}
    cells::Vector{C}
    nodes::Vector{N}
end

@testset "create_chunks with non-Grid AbstractGrid (BUG-013)" begin
    grid = generate_grid(Quadrilateral, (4, 4))
    dummygrid = DummyGrid(grid.cells, grid.nodes)
    cellset = collect(1:getncells(grid))

    # Automatic coloring path (colors_or_chunks = nothing)
    chunks_grid = FerriteAssembly.create_chunks(grid, cellset, nothing)
    chunks_dummy = FerriteAssembly.create_chunks(dummygrid, cellset, nothing)
    @test Set(Iterators.flatten(Iterators.flatten(chunks_dummy))) == Set(cellset)
    @test [sort!(collect(Iterators.flatten(c))) for c in chunks_dummy] == [sort!(collect(Iterators.flatten(c))) for c in chunks_grid]

    # User-supplied chunks path (doesn't use the grid argument beyond dispatch)
    chunks = [[cellset[1:8], cellset[9:16]], [cellset[17:end]]]
    @test FerriteAssembly.create_chunks(dummygrid, cellset, chunks) == chunks
end

@testset "work! with empty custom chunks (PR89)" begin
    # End-to-end reproduction: a domain whose user-supplied chunks contain
    # empty sub-chunks must still visit every cell during threaded `work!`.
    grid = generate_grid(Quadrilateral, (2, 1))
    ip = Lagrange{RefQuadrilateral, 1}()
    dh = close!(add!(DofHandler(grid), :u, ip))
    cv = CellValues(QuadratureRule{RefQuadrilateral}(2), ip)
    for num_tasks in (1, 2)
        leading_empty = [Int[] for _ in 1:num_tasks]
        chunks = [[leading_empty..., [1], Int[], [2], Int[]]]
        db = setup_domainbuffer(
            DomainSpec(dh, nothing, cv; chunks=chunks);
            threading=true, num_tasks=num_tasks)
        ig = SimpleIntegrator(Returns(1.0), 0.0)
        work!(ig, db)
        @test ig.val ≈ 4.0 # total area of the two-cell grid
    end
end