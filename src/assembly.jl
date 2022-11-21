
"""
    create_threaded_assemblers(K, r; nthreads=Threads.nthreads())

Convenience function for creating a `nthreads` long vector with the 
output of `start_assemble` as elements
"""
create_threaded_assemblers(K, r; nthreads=Threads.nthreads()) = [start_assemble(K,r) for _ in 1:nthreads]

""" 
    doassemble!(
        assembler::Ferrite.AbstractSparseAssembler, cellbuffer::AbstractCellBuffer, 
        s::AbstractVector, dh::DofHandler, 
        a=nothing, aold=nothing, Δt=nothing
        )

Sequential assembly of cells with the `dh::DofHandler`.
* `assembler` is obtained from `Ferrite.jl`'s `start_assemble(K,r)` function
* `cellbuffer` contains buffers for the specific cell. 
  See  [`CellBuffer`](@ref) for more info. 
* `s` is a collection (vector, dict, etc.) of state variables, where indexing 
  by `cellnr` gives the state variables for that cell. 
* `a` and `aold` are the current and old unknowns (can be set to `nothing` if not used)
* `Δt` is the time increment passed into each element routine
"""
function doassemble!(
    assembler::Ferrite.AbstractSparseAssembler, cellbuffer::AbstractCellBuffer, 
    states, dh::DofHandler, a=nothing, aold=nothing, Δt=nothing; cellset=1:Ferrite.getncells(dh)
    )
    for cellnr in cellset
        try
            assemble_cell!(assembler, cellbuffer, dh, cellnr, a, aold, states[cellnr], Δt)
        catch e
            rethrow(e)
        end
        
    end
end

""" 
    doassemble!(
        assemblers::Vector{<:Ferrite.AbstractSparseAssembler},
        cellbuffers::Vector{<:AbstractCellBuffer}, states, 
        dh::DofHandler, colored_sets::Vector{Vector{Int}}, 
        a=nothing, aold=nothing, Δt=nothing
        )

Threaded assembly of cells with the `dh::DofHandler`.
* `assemblers` are assemblers for each thread, which can be obtained
  with the [`create_threaded_assemblers`](@ref) function. 
* `cellbuffers` contains buffers for the specific cell,
  for each thread. This can be created by [`create_threaded_CellBuffers`](@ref). 
  See also [`CellBuffer`](@ref) for more info.
* `states`, `a`, `aold`, and `Δt` are the same as for the 
  sequential `doassemble!`
* `colored_sets` are cellsets for each color
"""
function doassemble!(
    assemblers::Vector{<:Ferrite.AbstractSparseAssembler},
    cellbuffers::Vector{<:AbstractCellBuffer}, states, 
    dh::DofHandler, colored_sets::Vector{Vector{Int}}, 
    a=nothing, aold=nothing, Δt=nothing; cellset=nothing
    )
    if length(assemblers) != length(cellbuffers) != Threads.nthreads()
        throw(DimensionMismatch("assemblers and cellbuffers must have as many elements as there are threads"))
    end
    for colorset in colored_sets
        Threads.@threads :static for cellnr in intersect_nothing(colorset, cellset)
            id = Threads.threadid()
            assemble_cell!(assemblers[id], cellbuffers[id], dh, cellnr, a, aold, states[cellnr], Δt)
        end
    end
end

""" 
    doassemble!(
        assembler::Ferrite.AbstractSparseAssembler, 
        cellbuffers::Tuple, states::Tuple, 
        dh::MixedDofHandler, 
        a=nothing, aold=nothing, Δt=nothing
        )

Sequential assembly of cells with the `dh::MixedDofHandler`.
* `cellbuffers` contains buffers for the specific cell in each `FieldHandler`
  in `dh.fieldhandlers` See  [`CellBuffer`](@ref) for more info. 
* `states` is a tuple which contains a collection of state variables, 
  one vector for each `FieldHandler` in `dh.fieldhandlers`. 
  Each element in the collection, i.e. `states[cellnr]`, contains the 
  state variables for one `Cell` with global number `cellnr`. 
  Unless all cells have the same type of the state, it might make sense to use a 
  Dict{Int,State} where the key refers to the global number. 
* `assembler`, `a`, `aold`, and `Δt` are the same as for the `DofHandler` case. 
"""
function doassemble!(
    assembler::Ferrite.AbstractSparseAssembler, 
    cellbuffers::Tuple, states::Tuple, 
    dh::MixedDofHandler, 
    a=nothing, aold=nothing, Δt=nothing; kwargs...
    )
    for (fh, cellbuffer, state) in zip(dh.fieldhandlers, cellbuffers, states)
        inner_doassemble!(assembler, cellbuffer, state, dh, fh, a, aold, Δt; kwargs...)
    end
end

""" 
    doassemble!(
        assemblers::Vector{<:Ferrite.AbstractSparseAssembler},
        cellbuffers::Tuple, 
        states::Tuple, 
        dh::MixedDofHandler, colored_sets::Vector{Vector{Int}}, 
        a::AbstractVector, aold::AbstractVector, Δt::Number
        )

Threaded assembly of cells with the `dh::MixedDofHandler`.
* `assemblers` and `colored_sets` are the same as for the threaded `DofHandler` case.
* `states` are the same as for the sequential `MixedDofHandler` case.
* `cellbuffers` contains vectors `Vector{AbstractCellBuffer}` for the cell type in 
  each `FieldHandler` in `dh.fieldhandlers`. The vector element corresponds to each 
  thread. This can be created by [`create_threaded_CellBuffers`](@ref). 
  See also [`CellBuffer`](@ref) for more info.
* `a`, `aold`, and `Δt` are the same as for the `DofHandler` case.
"""
function doassemble!(
    assemblers::Vector{<:Ferrite.AbstractSparseAssembler},
    cellbuffers::Tuple, 
    states::Tuple, 
    dh::MixedDofHandler, colored_sets::Vector{Vector{Int}}, 
    a=nothing, aold=nothing, Δt=nothing; kwargs...
    )
    for (fh, cellbuffer, state) in zip(dh.fieldhandlers, cellbuffers, states)
        inner_doassemble!(assemblers, cellbuffer, state, dh, fh, colored_sets, a, aold, Δt; kwargs...)
    end
end

function doassemble!(
    assemblers::Union{Vector{<:Ferrite.AbstractSparseAssembler}, Ferrite.AbstractSparseAssembler},
    cellbuffers::Dict{String},
    states::Dict{String},
    dh::Ferrite.AbstractDofHandler,
    args...)
    for (key,cellbuffer) in cellbuffers
        doassemble!(assemblers, cellbuffer, states[key], dh, args...; cellset=getcellset(dh, key))
    end
end

"""
    inner_doassemble!(
        assembler, cellbuffer::AbstractCellBuffer, states, 
        dh::MixedDofHandler, fh::FieldHandler, a, aold, Δt
        )

Sequential assembly of cells corresponding to the given `fh` 
from `dh.fieldhandlers`. 
Internal function that is called from the sequential version 
of `doassemble!` for the `MixedDofHandler`
"""
function inner_doassemble!(
    assembler, cellbuffer::AbstractCellBuffer, 
    states, dh::MixedDofHandler, fh::FieldHandler, 
    a, aold, Δt; cellset=nothing
    )
    for cellnr in intersect_nothing(fh.cellset, cellset)
        assemble_cell!(assembler, cellbuffer, dh, fh, cellnr, a, aold, states[cellnr], Δt)
    end
end

"""
    inner_doassemble!(
        assemblers::Vector{<:Ferrite.AbstractSparseAssembler},
        cellbuffers::Vector{<:AbstractCellBuffer}, states, 
        dh::MixedDofHandler, fh::FieldHandler, 
        colored_sets::Vector{Vector{Int}}, 
        a, aold, Δt)

Parallel assembly of cells corresponding to the given `fh` from 
`dh.fieldhandlers`.
Internal function that is called from the parallel version of 
`doassemble!` for the `MixedDofHandler`
"""
function inner_doassemble!(
    assemblers::Vector{<:Ferrite.AbstractSparseAssembler},
    cellbuffers::Vector{<:AbstractCellBuffer}, states, 
    dh::MixedDofHandler, fh::FieldHandler, colored_sets::Vector{Vector{Int}}, 
    a, aold, Δt; cellset=nothing)

    if length(assemblers) != length(cellbuffers) != Threads.nthreads()
        throw(DimensionMismatch("assemblers and cellbuffers must have as many elements as there are threads"))
    end
    for colorset in colored_sets
        Threads.@threads :static for cellnr in intersect_nothing(intersect(colorset, fh.cellset), cellset)
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
    assemble_cell_reinited!(assembler, cellbuffer, dh_fh::Union{DofHandler,FieldHandler}, state, Δt)

Internal function that assembles the cell described by the reinitialized `cellbuffer`. This function is called 
in all cases: Parallel or sequential and `DofHandler` or `MixedDofHandler`
"""
function assemble_cell_reinited!(assembler, cellbuffer, dh_fh::Union{DofHandler,FieldHandler}, state, Δt)
    Ke = get_Ke(cellbuffer)
    re = get_re(cellbuffer)
    ae = get_ae(cellbuffer)
    material = get_material(cellbuffer)
    cellvalues = get_cellvalues(cellbuffer)
    element_routine!(Ke, re, state, ae, material, cellvalues, dh_fh, Δt, cellbuffer)
    dofs = celldofs(cellbuffer)
    assemble!(assembler, dofs, Ke, re)
end