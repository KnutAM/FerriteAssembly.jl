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