# Inspired by Base.Channel, but data is not removed, 
# only the index is incremented. 
mutable struct ChunkIterator{T}
    const lock::ReentrantLock
    const chunks::Vector{Vector{T}}
    index::Int
end
ChunkIterator(chunks::Vector{<:Vector}) = ChunkIterator(ReentrantLock(), chunks, 0)
reset!(ci::ChunkIterator) = (ci.index=0)

# Base.AbstractLock interface
Base.lock(ci::ChunkIterator) = lock(ci.lock)
Base.lock(f, ci::ChunkIterator) = lock(f, ci.lock)
Base.unlock(ci::ChunkIterator) = unlock(ci.lock)
Base.trylock(ci::ChunkIterator) = trylock(ci.lock)

# Iterator interface
function Base.iterate(ci::ChunkIterator, state=nothing)
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
Base.length(ci::ChunkIterator) = length(ci.chunks)
Base.eltype(::ChunkIterator{T}) where T = Vector{T}

# Non-iterator implementation
function get_chunk(ci::ChunkIterator{T}) where T
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