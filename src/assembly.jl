setup_assemblers(::Val{false}, K, r; fillzero=true) = start_assemble(K, r; fillzero=fillzero)
function setup_assemblers(::Val{true}, K, r; fillzero=true)
    assemblers = [start_assemble(K, r; fillzero=fillzero)]
    for _ in 2:Threads.nthreads()
        push!(assemblers, start_assemble(K, r; fillzero=false))
    end
    return assemblers
end

function doassemble!(K::AbstractMatrix, r::AbstractVector, new_states::Dict{Int}, old_states::Dict{Int}, buffer::DomainBuffer; kwargs...)
    doassemble!(K, r, Dict("noname"=>new_states), Dict("noname"=>old_states), Dict("noname"=>buffer); kwargs...)
end
function doassemble!(r::AbstractVector{<:Number}, new_states::Dict{Int}, old_states::Dict{Int}, buffer::DomainBuffer; kwargs...)
    doassemble!(r, Dict("noname"=>new_states), Dict("noname"=>old_states), Dict("noname"=>buffer); kwargs...)
end

function doassemble!(K::AbstractMatrix, r::AbstractVector, new_states::Dict{String}, old_states::Dict{String}, buffers::Dict{String}; fillzero=true, kwargs...)
    threaded = Val(iscolored(buffers))
    assemblers = setup_assemblers(threaded, K, r; fillzero=fillzero)
    _doassemble!(threaded, assemblers, new_states, old_states, buffers; kwargs...)
end
function doassemble!(r::AbstractVector{<:Number}, new_states::Dict{String}, old_states::Dict{String}, buffers::Dict{String}; fillzero=true, kwargs...)
    threaded = Val(iscolored(buffers))
    fillzero && fill!(r, zero(eltype(r)))
    _doassemble!(threaded, r, new_states, old_states, buffers; kwargs...)
end

function _doassemble!(threaded::Val, assembler_or_r, new_states, old_states, buffers; a=nothing, aold=nothing, Δt=NaN)
    buffer = first(values(buffers))
    reset_scaling!(buffer.scaling)
    reset_scaling!.(buffer.scalings)
    for key in keys(new_states)
        buffer = buffers[key]
        update_time!(buffer.cellbuffer, Δt)
        assemble_domain!(threaded, assembler_or_r, new_states[key], old_states[key], buffer, a, aold)
    end
    for scaling in buffer.scalings
        add_to_scaling!(buffer.scaling, scaling)
    end
end

function assemble_domain!(#=threaded=#::Val{false}, assembler_or_r, new_states, old_states, buffer, a, aold)
    for cellnr in buffer.cellset
        assemble_cell!(assembler_or_r, buffer.cellbuffer, buffer.sdh, cellnr, a, aold, new_states, old_states, buffer.scaling)
    end
end

function assemble_domain!(#=threaded=#::Val{true}, assemblers::Vector{<:Ferrite.AbstractSparseAssembler}, new_states, old_states, buffer, a, aold)
    cellbuffers = buffer.cellbuffer
    scalings = buffer.scalings
    for cellset in buffer.cellset
        Threads.@threads :static for cellnr in cellset
            id = Threads.threadid()
            assemble_cell!(assemblers[id], cellbuffers[id], buffer.sdh, cellnr, a, aold, new_states, old_states, scalings[id])
        end
    end
end

function assemble_domain!(#=threaded=#::Val{true}, r::Vector{<:Number}, new_states, old_states, buffer, a, aold)
    cellbuffers = buffer.cellbuffer
    scalings = buffer.scalings
    for cellset in buffer.cellset
        Threads.@threads :static for cellnr in cellset
            id = Threads.threadid()
            assemble_cell!(r, cellbuffers[id], buffer.sdh, cellnr, a, aold, new_states, old_states, scalings[id])
        end
    end
end

"""
    assemble_cell!(assembler, cellbuffer, dh::DofHandler, cellnr, a, aold, state, Δt)
    assemble_cell!(assembler, cellbuffer, dh::MixedDofHandler, fh::FieldHandler, cellnr, a, aold, state, Δt)

Internal function to that reinitializes the `cellbuffer` and calls [`assemble_cell_reinited!`](@ref).
"""
function assemble_cell!(assembler_or_r, cellbuffer, sdh::SubDofHandler, cellnr, a, aold, new_states, old_states, scaling)
    reinit!(cellbuffer, sdh.dh, cellnr, a, aold, old_states)
    cell_state = fast_getindex(new_states, cellnr)
    assemble_cell_reinited!(assembler_or_r, cellbuffer, cell_state, scaling)
end

"""
    assemble_cell_reinited!(assembler, cellbuffer, state, scaling)

Internal function that assembles the cell described by the reinitialized `cellbuffer`.
"""
function assemble_cell_reinited!(assembler::Ferrite.AbstractSparseAssembler, cellbuffer, state, scaling)
    Ke = get_Ke(cellbuffer)
    re = get_re(cellbuffer)
    ae = get_ae(cellbuffer)
    material = get_material(cellbuffer)
    cellvalues = get_cellvalues(cellbuffer)
    element_routine!(Ke, re, state, ae, material, cellvalues, cellbuffer)
    update_scaling!(scaling, re, cellbuffer)
    dofs = celldofs(cellbuffer)
    assemble!(assembler, dofs, Ke, re)
end

"""
    assemble_cell_reinited!(assembler, cellbuffer, state, scaling)

Internal function that assembles the residual for the cell described by the reinitialized `cellbuffer`.
"""
function assemble_cell_reinited!(r::Vector{<:Number}, cellbuffer, state, scaling)
    re = get_re(cellbuffer)
    ae = get_ae(cellbuffer)
    material = get_material(cellbuffer)
    cellvalues = get_cellvalues(cellbuffer)
    element_residual!(re, state, ae, material, cellvalues, cellbuffer)
    update_scaling!(scaling, re, cellbuffer)
    dofs = celldofs(cellbuffer)
    # assemble!(r, dofs, re) would be nice. 
    for (i, d) in pairs(dofs)
        r[d] += re[i]
    end
end
