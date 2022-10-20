
"""
    create_threaded_assemblers(K, r; nthreads=Threads.nthreads())

Convenience function for creating a `nthreads` long vector with the 
output of `start_assemble` as elements
"""
create_threaded_assemblers(K, r; nthreads=Threads.nthreads()) = [start_assemble(K,r) for _ in 1:nthreads]

""" 
    doassemble!(
        assembler::Ferrite.AbstractSparseAssembler, cellbuffer::CellBuffer, 
        s::AbstractVector, dh::DofHandler, 
        a::AbstractVector, aold::AbstractVector, Δt::Number
        )

Sequential assembly of cells with the `dh::DofHandler`.
* `assembler` is obtained from `Ferrite.jl`'s `start_assemble(K,r)` function
* `cellbuffer` contains buffers for the specific cell. 
  See  [`CellBuffer`](@ref) for more info. 
* `s` is a vector of state variables, where each
  element contains state variables for each `cellnr`. 
* `a` and `aold` are the current and old unknowns. 
* `Δt` is the time increment passed into each element routine
"""
function doassemble!(
    assembler::Ferrite.AbstractSparseAssembler, cellbuffer::CellBuffer, 
    s::AbstractVector, dh::DofHandler, 
    a::AbstractVector, aold::AbstractVector, Δt::Number
    )
    for cellnr in 1:Ferrite.getncells(dh)
        assemble_cell!(assembler, cellbuffer, dh, cellnr, a, aold, s[cellnr], Δt)
    end
end

""" 
    doassemble!(
        assemblers::Vector{<:Ferrite.AbstractSparseAssembler},
        cellbuffers::Vector{<:CellBuffer}, 
        s::AbstractVector, 
        colored_sets::Vector{Vector{Int}}, dh::DofHandler, 
        a::AbstractVector, aold::AbstractVector, Δt::Number
        )

Threaded assembly of cells with the `dh::DofHandler`.
* `assemblers` are assemblers for each thread, which can be obtained
  with the [`create_threaded_assemblers`](@ref) function. 
* `cellbuffers` contains buffers for the specific cell,
  for each thread. This can be created by [`create_threaded_CellBuffers`](@ref). 
  See also [`CellBuffer`](@ref) for more info.
* `s`, `a`, `aold`, and `Δt` are the same as for the 
  sequential `doassemble!`
* `colored_sets` are cellsets for each color
"""
function doassemble!(
    assemblers::Vector{<:Ferrite.AbstractSparseAssembler},
    cellbuffers::Vector{<:CellBuffer}, 
    s::AbstractVector, 
    colored_sets::Vector{Vector{Int}}, dh::DofHandler, 
    a::AbstractVector, aold::AbstractVector, Δt::Number
    )
    if length(assemblers) != length(cellbuffers) != Threads.nthreads()
        throw(DimensionMismatch("assemblers and cellbuffers must have as many elements as there are threads"))
    end
    for cellset in colored_sets
        Threads.@threads for cellnr in cellset
            id = Threads.threadid()
            assemble_cell!(assemblers[id], cellbuffers[id], dh, cellnr, a, aold, s[cellnr], Δt)
        end
    end
end

""" 
    doassemble!(
        assembler::Ferrite.AbstractSparseAssembler, 
        cellbuffers::Tuple, states::Tuple, 
        dh::MixedDofHandler, 
        a::AbstractVector, aold::AbstractVector, Δt::Number
        )

Sequential assembly of cells with the `dh::MixedDofHandler`.
* `cellbuffers` contains buffers for the specific cell in each `FieldHandler`
  in `dh.fieldhandlers` See  [`CellBuffer`](@ref) for more info. 
* `states` is a tuple which contains vectors of state variables, 
  one vector for each `FieldHandler` in `dh.fieldhandlers`. 
  Each vector element contains the state variables for one `Cell`. 
  Note that the vector index corresponds to the cellnr in the grid, 
  and not in the cellset of the `FieldHandler` 
  (as this is a `Set` and the order is not guaranteed). Unless all cells
  have the same type of the state, it might make sense to use a sparse vector. 
* `assembler`, `a`, `aold`, and `Δt` are the same as for the `DofHandler` case. 
"""
function doassemble!(
    assembler::Ferrite.AbstractSparseAssembler, 
    cellbuffers::Tuple, states::Tuple, 
    dh::MixedDofHandler, 
    a::AbstractVector, aold::AbstractVector, Δt::Number
    )
    for (fh, cellbuffer, state) in zip(dh.fieldhandlers, cellbuffers, states)
        inner_doassemble!(assembler, cellbuffer, state, dh, fh, a, aold, Δt)
    end
end

""" 
    doassemble!(
        assemblers::Vector{<:Ferrite.AbstractSparseAssembler},
        cellbuffers::Tuple, 
        states::Tuple, 
        colored_sets::Vector{Vector{Int}}, dh::MixedDofHandler, 
        a::AbstractVector, aold::AbstractVector, Δt::Number
        )

Threaded assembly of cells with the `dh::MixedDofHandler`.
* `assemblers` and `colored_sets` are the same as for the threaded `DofHandler` case.
* `states` are the same as for the sequential `MixedDofHandler` case.
* `cellbuffers` contains vectors `Vector{CellBuffer}` for the cell type in 
  each `FieldHandler` in `dh.fieldhandlers`. The vector element corresponds to each 
  thread. This can be created by [`create_threaded_CellBuffers`](@ref). 
  See also [`CellBuffer`](@ref) for more info.
* `a`, `aold`, and `Δt` are the same as for the `DofHandler` case.
"""
function doassemble!(
    assemblers::Vector{<:Ferrite.AbstractSparseAssembler},
    cellbuffers::Tuple, 
    states::Tuple, 
    colored_sets::Vector{Vector{Int}}, dh::MixedDofHandler, 
    a::AbstractVector, aold::AbstractVector, Δt::Number
    )
    for (fh, cellbuffer, state) in zip(dh.fieldhandlers, cellbuffers, states)
        inner_doassemble!(assemblers, cellbuffer, state, colored_sets, dh, fh, a, aold, Δt)
    end
end

"""
    inner_doassemble!(
        assembler, cellbuffer::CellBuffer, state, 
        dh::MixedDofHandler, fh::FieldHandler, a, aold, Δt
        )

Sequential assembly of cells corresponding to the given `fh` 
from `dh.fieldhandlers`. 
Internal function that is called from the sequential version 
of `doassemble!` for the `MixedDofHandler`
"""
function inner_doassemble!(
    assembler, cellbuffer::CellBuffer, 
    state, dh::MixedDofHandler, fh::FieldHandler, 
    a, aold, Δt
    )
    for cellnr in fh.cellset
        assemble_cell!(assembler, cellbuffer, dh, fh, cellnr, a, aold, state[cellnr], Δt)
    end
end

"""
    inner_doassemble!(
        assemblers::Vector{<:Ferrite.AbstractSparseAssembler},
        cellbuffers::Vector{<:CellBuffer}, 
        states::AbstractVector, colored_sets::Vector{Vector{Int}}, 
        dh::MixedDofHandler, fh::FieldHandler, a, aold, Δt)

Parallel assembly of cells corresponding to the given `fh` from 
`dh.fieldhandlers`.
Internal function that is called from the parallel version of 
`doassemble!` for the `MixedDofHandler`
"""
function inner_doassemble!(
    assemblers::Vector{<:Ferrite.AbstractSparseAssembler},
    cellbuffers::Vector{<:CellBuffer}, 
    states::AbstractVector, 
    colored_sets::Vector{Vector{Int}}, dh::MixedDofHandler, fh::FieldHandler, 
    a, aold, Δt)

    if length(assemblers) != length(cellbuffers) != Threads.nthreads()
        throw(DimensionMismatch("assemblers and cellbuffers must have as many elements as there are threads"))
    end
    for cellset in colored_sets
        Threads.@threads for cellnr in intersect(cellset, fh.cellset)
            id = Threads.threadid()
            assemble_cell!(assemblers[id], cellbuffers[id], dh, fh, cellnr, a, aold, states[cellnr], Δt)
        end
    end
end

"""
    assemble_cell!(assembler, cellbuffer, dh::DofHandler, cellnr, a, aold, state, Δt)
    assemble_cell!(assembler, cellbuffer, dh::MixedDofHandler, fh::FieldHandler, cellnr, a, aold, state, Δt)

Internal function to that reinitializes the `cellbuffer` and calls [`assemble_cell_reinited!`](@ref).
"""
function assemble_cell!(assembler, cellbuffer, dh::DofHandler, cellnr, a, aold, state, Δt)
    reinit!(cellbuffer, dh, cellnr, a, aold)
    assemble_cell_reinited!(assembler, cellbuffer, dh, state, Δt)
end

function assemble_cell!(assembler, cellbuffer, dh::MixedDofHandler, fh::FieldHandler, cellnr, a, aold, state, Δt)
    reinit!(cellbuffer, dh, cellnr, a, aold)
    assemble_cell_reinited!(assembler, cellbuffer, fh, state, Δt)
end 


"""
    assemble_cell_reinited!(assembler, cellbuffer::CellBuffer, dh_fh::Union{DofHandler,FieldHandler}, state, Δt)

Internal function that assembles the cell described by the reinitialized `cellbuffer`. This function is called 
in all cases: Parallel or sequential and `DofHandler` or `MixedDofHandler`
"""
function assemble_cell_reinited!(assembler, cellbuffer::CellBuffer, dh_fh::Union{DofHandler,FieldHandler}, state, Δt)
    (;Ke, re, ae, ae_old, material, cellvalues, cache, dofs) = cellbuffer

    unscale_primary!.((ae, ae_old), (material,), (dh_fh,))
    element_routine!(
        Ke, re, ae, ae_old, state, material, 
        cellvalues, dh_fh, Δt, cache
        )
    scale_residual!(Ke, re, material, dh_fh)

    assemble!(assembler, dofs, Ke, re)
end

