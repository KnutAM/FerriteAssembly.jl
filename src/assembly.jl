# Assemble single domain
function doassemble!(assembler, new_states, buffer::DomainBuffer; a=nothing, aold=nothing, old_states=nothing, Δt=NaN)
    assemble_domain!(assembler, new_states, buffer, a, aold, old_states, Δt)
end
function doassemble!(assembler, new_states, buffer::ThreadedDomainBuffer; a=nothing, aold=nothing, old_states=nothing, Δt=NaN)
    threaded_assembler = TaskLocals(assembler)
    assemble_domain!(threaded_assembler, new_states, buffer, a, aold, old_states, Δt)
end
function doassemble!(assembler::TaskLocals, new_states, buffer::ThreadedDomainBuffer; a=nothing, aold=nothing, old_states=nothing, Δt=NaN)
    assemble_domain!(assembler, new_states, buffer, a, aold, old_states, Δt)
end

# Assemble multiple domains
function doassemble!(assembler, new_states, buffers::Dict{String,<:DomainBuffer}; a=nothing, aold=nothing, old_states=nothing, Δt=NaN)
    for key in keys(new_states)
        domain_old_states = isnothing(old_states) ? nothing : old_states[key]
        assemble_domain!(assembler, new_states[key], buffers[key], a, aold, domain_old_states, Δt)
    end
end
function doassemble!(assembler, new_states, buffers::Dict{String,<:ThreadedDomainBuffer}; kwargs...)
    threaded_assembler = TaskLocals(assembler)
    doassemble!(threaded_assembler, new_states, buffers; kwargs...)
end
function doassemble!(assembler::TaskLocals, new_states, buffers::Dict{String,<:ThreadedDomainBuffer}; a=nothing, aold=nothing, old_states=nothing, Δt=NaN)
    for key in keys(new_states)
        domain_old_states = isnothing(old_states) ? nothing : old_states[key]
        assemble_domain!(assembler, new_states[key], buffers[key], a, aold, domain_old_states, Δt)
    end
end

function assemble_domain!(assembler, new_states, buffer::DomainBuffer, a, aold, old_states, Δt)
    cellbuffer = buffer.cellbuffer
    update_time!(cellbuffer, Δt)
    for cellnr in buffer.cellset
        assemble_cell!(assembler, new_states, cellbuffer, buffer.sdh, cellnr, a, aold, old_states)
    end
end

function assemble_domain!(assembler::TaskLocals, new_states, buffer::ThreadedDomainBuffer, a, aold, old_states, Δt)
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

"""
    assemble_cell_reinited!(assembler::Ferrite.AbstractSparseAssembler, state, cellbuffer)

Internal function that assembles local element stiffness, `Ke`, and the local residual, `re`, 
for the cell described by the reinitialized `cellbuffer` into the global matrix and vector in `assembler`
"""
function assemble_cell_reinited!(assembler::Ferrite.AbstractSparseAssembler, state, cellbuffer)
    Ke = get_Ke(cellbuffer)
    re = get_re(cellbuffer)
    ae = get_ae(cellbuffer)
    material = get_material(cellbuffer)
    cellvalues = get_cellvalues(cellbuffer)
    element_routine!(Ke, re, state, ae, material, cellvalues, cellbuffer)
    assemble!(assembler, celldofs(cellbuffer), Ke, re)
end

"""
    assemble_cell_reinited!(assembler::KeReAssembler, state, cellbuffer)

Internal function that assembles local element stiffness, `Ke`, and the local residual, `re`, 
for the cell described by the reinitialized `cellbuffer` into the global matrix and vector in 
`assembler`. In addition, `KeReAssembler` supports extra options, in particular
- Local application of constraints Ferrite.apply_assemble!
- Residual scaling factor, e.g. `ElementResidualScaling`
"""
function assemble_cell_reinited!(assembler::KeReAssembler, state, cellbuffer)
    Ke = get_Ke(cellbuffer)
    re = get_re(cellbuffer)
    ae = get_ae(cellbuffer)
    material = get_material(cellbuffer)
    cellvalues = get_cellvalues(cellbuffer)
    element_routine!(Ke, re, state, ae, material, cellvalues, cellbuffer)
    assemble!(assembler, cellbuffer)
end

"""
    assemble_cell_reinited!(assembler::ReAssembler, cellbuffer, state, scaling)

Internal function that assembles the residual for the cell described by the reinitialized `cellbuffer`.
"""
function assemble_cell_reinited!(assembler::ReAssembler, state, cellbuffer)
    re = get_re(cellbuffer)
    ae = get_ae(cellbuffer)
    material = get_material(cellbuffer)
    cellvalues = get_cellvalues(cellbuffer)
    element_residual!(re, state, ae, material, cellvalues, cellbuffer)
    
    assemble!(assembler, cellbuffer)
end
