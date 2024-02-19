# Inspired by Base.Channel, but data is not removed, 
# only the index is incremented. 
"""
    TaskChunks(chunks)

This is **strictly** an internal implementation detail for now.
It is used to provide a variant of `Base.Channel`, but 
where data is not removed, the index is just incremented.
"""
mutable struct TaskChunks{T}
    const lock::ReentrantLock
    const chunks::Vector{Vector{T}}
    index::Int
end
TaskChunks(chunks::Vector{<:Vector}) = TaskChunks(ReentrantLock(), chunks, 0)

# Base.AbstractLock interface
Base.lock(ci::TaskChunks) = lock(ci.lock)
Base.unlock(ci::TaskChunks) = unlock(ci.lock)
# These are not required. 
#Base.lock(f, ci::TaskChunks) = lock(f, ci.lock)
#Base.trylock(ci::TaskChunks) = trylock(ci.lock)

# Non-iterator implementation
function get_chunk(ci::TaskChunks{T}) where T
    lock(ci)
    try
        if length(ci.chunks) <= ci.index
            return T[]
        else
            ci.index += 1
            return @inbounds ci.chunks[ci.index]
        end
    finally
        unlock(ci)
    end
end

# Create chunks 
function split_in_chunks(set::Vector{T}; num_tasks = Threads.nthreads()) where T #::Vector{Vector{T}}
    # Split `set` into chunks of indices to be assembled
    # Each task takes one chunk if indices, assemble them, and asks for a new chunk
    # The "standard" threading behavior is obtained by making one chunk per thread,
    # this would be ideal if all indices take equal time to assemble. 
    # The opposite would be one index in each chunk, which ensures very good load distribution.
    # However, overhead is large as asking for new chunks require a lock to prevent race conditions. 

    # This function tries to make a reasonable default chunking and 
    # accounts for the number of cells and threads.
    # It is called when setting up the ThreadedDomainBuffer.
    # Some checks are made here in addition to unit tests to ensure correct distribution.
    # Errors in that would otherwise cause wrong values silently during assembly.
    num_items = length(set)
    if num_items < num_tasks
        chunk_size = 0
    elseif num_tasks == 1
        chunk_size = num_items
    else
        max_chunk_size = 100
        # Not carefully evaluated heuristic, but gives somewhat reasonable values by looking at them
        chunk_size = max(min(max_chunk_size, round(Int,sqrt(num_items/(2*num_tasks)))), 1)
    end
    num_chunks = num_items ÷ (chunk_size+1) + 1
    num_missing = num_items - num_chunks*chunk_size # How many chunks that will not have the full number 
    # Check to be sure
    num_items == (chunk_size*num_chunks + num_missing) || error("This is a bug")
    num_missing <= num_chunks || error("This is a bug")
    chunks = Vector{T}[]
    i1 = 1
    for _ in 1:num_chunks
        Δi = num_missing>0 ? 1 : 0 # To distribute the lower numbers
        i2 = i1 + (chunk_size-1) + Δi
        push!(chunks, set[i1:i2])
        i1 = i2+1
        num_missing = max(num_missing-1, 0)
    end
    # Check to be sure
    sum(length, chunks)==num_items || error("This is a bug")
    return chunks
end

# If chunks already given, check that they match the intersected set
function create_chunks(::Grid, intersected_set::Vector{I}, chunks::Vector{Vector{Vector{I}}}) where I
    chunk_set = sizehint!(Set{I}(), length(intersected_set))
    full_set = Set(intersected_set)
    for chunk_vector in chunks
        for chunk in chunk_vector
            for cellnr in chunk
                if cellnr ∈ chunk_set
                    error("cell ", cellnr, " in the chunks occurs at least twice")
                elseif cellnr ∉ full_set
                    error("cell ", cellnr, " in chunks are not in the set")
                end
                push!(chunk_set, cellnr)
            end
        end
    end
    full_set==chunk_set || error("Not all items in set where included in the chunks")
    return chunks
end

# Colors given with same type as set 
function create_chunks(g::Grid, intersected_set::Vector{I}, colors::Vector{Vector{I}}) where I
    colors_intersect = map(sort! ∘ collect ∘ Base.Fix1(intersect, intersected_set), colors)
    chunks = [split_in_chunks(set) for set in colors_intersect]
    return create_chunks(g, intersected_set, chunks)
end
# Colors given with different type as set
function create_chunks(g::Grid, intersected_set::Vector{FaceIndex}, colors::Vector{Vector{Int}})
    cellset = first.(intersected_set)
    colors_intersect = map(sort! ∘ collect ∘ Base.Fix1(intersect, cellset), colors)
    chunks = [convert_chunk(split_in_chunks(set), intersected_set) for set in colors_intersect]
    return create_chunks(g, intersected_set, chunks) # Pass through to check correctness
end

# Colors not given
function create_chunks(grid::Grid, intersected_set::Vector, ::Nothing)
    makecellset(v::Vector{Int}) = Set(v)
    makecellset(v::Vector) = Set(first.(v))
    return create_chunks(grid, intersected_set, create_coloring(grid, makecellset(intersected_set)))
end

# Convert chunks for cells to chunks for faces by adding all face::FaceIndex with the same cell in the same chunk
function convert_chunk(cellchunks::Vector{Vector{Int}}, set::Vector{FaceIndex})
    facechunks = Vector{FaceIndex}[]
    facechunk = FaceIndex[] # workspace
    for cellchunk in cellchunks
        facechunk = sizehint!(FaceIndex[], length(cellchunk))
        for cellnr in cellchunk
            for (cellnr_set, facenr) in set
                cellnr_set == cellnr && push!(facechunk, FaceIndex(cellnr, facenr))
            end
        end
        push!(facechunks, facechunk)
    end
    return facechunks
end