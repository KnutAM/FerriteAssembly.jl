function create_chunks(set::Vector{Int}; num_tasks = Threads.nthreads()) #::Vector{Vector{Int}}
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
    num_cells = length(set)
    if num_cells < num_tasks
        chunk_size = 0
    elseif num_tasks == 1
        chunk_size = num_cells
    else
        max_chunk_size = 100
        # Not carefully evaluated heuristic, but gives somewhat reasonable values by looking at them
        chunk_size = max(min(max_chunk_size, round(Int,sqrt(num_cells/(2*num_tasks)))), 1)
    end
    num_chunks = num_cells÷(chunk_size+1) + 1
    num_missing = num_cells - num_chunks*chunk_size # How many chunks that will not have the full number 
    # Check to be sure
    if num_cells != (chunk_size*num_chunks + num_missing) || (num_missing > num_chunks)
        @show num_cells, num_tasks, chunk_size, num_chunks, num_missing
        error("This should not happen and is a bug")
    end
    chunks = Vector{Int}[]
    i1 = 1
    for _ in 1:num_chunks
        Δi = num_missing>0 ? 1 : 0 # To distribute the lower numbers
        i2 = i1 + (chunk_size-1) + Δi
        push!(chunks, set[i1:i2])
        i1 = i2+1
        num_missing = max(num_missing-1, 0)
    end
    # Check to be sure
    sum(length, chunks)==num_cells || error("This should not happen and is a bug")
    return chunks
end

# Inspired by Base.Channel, but data is not removed, 
# only the index is incremented. 
"""
    TaskChunks(chunks)

This is strictly an internal implementation detail for now.
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

# Iterator interface (not used for now)
#=
function Base.iterate(ci::TaskChunks, state=nothing)
    lock(ci)
    local retval
    try
        if length(ci.chunks) <= ci.index
            retval = nothing 
        else
            ci.index += 1
            retval = (@inbounds ci.chunks[ci.index], nothing)
        end
    finally
        unlock(ci)
    end
    return retval
end
Base.length(ci::TaskChunks) = length(ci.chunks)
Base.eltype(::TaskChunks{T}) where T = Vector{T}
=#