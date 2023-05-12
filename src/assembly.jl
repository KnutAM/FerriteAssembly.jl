"""
    doassemble!(assembler, new_states::Dict{Int}, buffer::AbstractDomainBuffer; a=nothing, aold=nothing, old_states=nothing, Δt=NaN)

Use `assembler` to assemble a single domain described by `buffer`, and update `new_states` if dictated by the assembler. 

    doassemble!(assembler, new_states::Dict{String}, buffers::Dict{String,AbstractDomainBuffer}; a=nothing, aold=nothing, old_states=nothing, Δt=NaN)

Use `assembler` to assemble all domains described by `buffers`, and update `new_states` if dictated by the assembler.

The keyword arguments work as follows:
- `a, aold`: If these are `nothing`, then `NaN` values are given to the element routine for the corresponding entry.
- `old_states`: If `nothing`, the old state in the cellbuffer will not be updated before calling the element routine 
- `Δt`: Directly the value returned by `get_time_increment(cellbuffer)`
"""
function doassemble!(assembler, new_states, buffer::DomainBuffer; kwargs...)
    assemble_domain!(assembler, new_states, buffer; kwargs...)
end
function doassemble!(assembler, new_states, buffer::ThreadedDomainBuffer; kwargs...)
    threaded_assembler = TaskLocals(assembler)
    assemble_domain!(threaded_assembler, new_states, buffer; kwargs...)
end

# Assemble multiple domains
function doassemble!(assembler, new_states, buffers::Dict{String,<:DomainBuffer}; old_states=nothing, kwargs...)
    for key in keys(new_states)
        domain_old_states = isnothing(old_states) ? nothing : fast_getindex(old_states, key)
        assemble_domain!(assembler, new_states[key], buffers[key]; old_states=domain_old_states, kwargs...)
    end
end
function doassemble!(assembler, new_states, buffers::Dict{String,<:ThreadedDomainBuffer}; kwargs...)
    threaded_assembler = TaskLocals(assembler)
    threaded_doassemble!(threaded_assembler, new_states, buffers; kwargs...)
end
function threaded_doassemble!(assembler::TaskLocals, new_states, buffers::Dict{String,<:ThreadedDomainBuffer}; old_states=nothing, kwargs...)
    for key in keys(new_states)
        domain_old_states = isnothing(old_states) ? nothing : fast_getindex(old_states, key)
        assemble_domain!(assembler, new_states[key], buffers[key]; old_states=domain_old_states, kwargs...)
    end
end

function assemble_domain!(assembler, new_states, buffer::DomainBuffer; a=nothing, aold=nothing, old_states=nothing, Δt=NaN)
    cellbuffer = buffer.cellbuffer
    update_time!(cellbuffer, Δt)
    for cellnr in buffer.cellset
        assemble_cell!(assembler, new_states, cellbuffer, buffer.sdh, cellnr, a, aold, old_states)
    end
end

function assemble_domain!(assembler::TaskLocals, new_states, buffer::ThreadedDomainBuffer; a=nothing, aold=nothing, old_states=nothing, Δt=NaN)
    cellbuffers = buffer.cellbuffers #::TaskLocals
    update_time!(get_base(cellbuffers), Δt)
    scatter!(cellbuffers)
    scatter!(assembler)
    for cellset in buffer.colors
        Threads.@threads :static for cellnr in cellset
            id = Threads.threadid()
            assemble_cell!(get_local(assembler, id), new_states, get_local(cellbuffers, id), buffer.sdh, cellnr, a, aold, old_states)
        end
    end
    gather!(cellbuffers)
    gather!(assembler)
end

"""
    assemble_cell!(assembler, new_states, cellbuffer, sdh::SubDofHandler, cellnr, a, aold, old_states)

Internal function to that reinitializes the `cellbuffer` and calls [`assemble_cell_reinited!`](@ref).
"""
function assemble_cell!(assembler, new_states, cellbuffer, sdh::SubDofHandler, cellnr, a, aold, old_states)
    reinit!(cellbuffer, sdh.dh, cellnr, a, aold, old_states)
    cell_state = fast_getindex(new_states, cellnr)
    assemble_cell_reinited!(assembler, cell_state, cellbuffer)
end