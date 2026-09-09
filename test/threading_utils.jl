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

@testset "get_chunk with empty chunks (BUG-002)" begin
    # A legitimately empty chunk (e.g. supplied by a user via `chunks=...`, or
    # produced by an unlucky coloring/splitting) must not be mistaken by a
    # worker for queue exhaustion: only running out of chunks should stop a
    # worker from asking for more. Every task must therefore be able to reach
    # chunks with real work located after leading/interspersed empty chunks.
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

@testset "work! with empty custom chunks (BUG-002)" begin
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