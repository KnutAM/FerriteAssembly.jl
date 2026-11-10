"""
    CoupledCellBuffer(primary::CellBuffer, partner_buffers::NamedTuple)

Wraps a reader's own `primary::CellBuffer` together with references to its coupling
partners' plain `CellBuffer`s (`partner_buffers`, returned by [`get_coupled_buffer`](@ref)).

`partner_buffers` are never themselves `CoupledCellBuffer`s or `AutoDiffCellBuffer`s:
coupling is direct and nonrecursive, so a partner's own coupling (if any) is not exposed.

Constructed once per task at [`CoupledSimulations`](@ref) setup time; never rebuilt per cell
or per `work!` call. Does *not* hold the partner `Simulation`s needed to reinitialize
`partner_buffers`: those live once on the [`CoupledSimulation`](@ref) passed into
[`reinit_buffer!`](@ref) (defined in `Coupling.jl`, once that type exists), rather than being
duplicated into every task-local copy of this buffer.
"""
struct CoupledCellBuffer{CB<:CellBuffer, PB<:NamedTuple} <: AbstractCellBufferWrapper
    primary::CB
    partner_buffers::PB
end

@inline get_wrapped(cb::CoupledCellBuffer) = cb.primary

get_coupled_buffers(cb::CoupledCellBuffer) = cb.partner_buffers

function _replace_material_with(cb::CoupledCellBuffer, new_material)
    new_primary = _replace_material_with(cb.primary, new_material)
    return CoupledCellBuffer(new_primary, cb.partner_buffers)
end
