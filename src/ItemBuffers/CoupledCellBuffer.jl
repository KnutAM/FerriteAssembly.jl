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
struct CoupledCellBuffer{CB<:CellBuffer, PB<:NamedTuple} <: AbstractCellBuffer
    primary::CB
    partner_buffers::PB
end

for op = (:get_Ke, :get_re, :get_ae, :get_material, :get_values, :get_time_increment,
        :get_aeold, :get_state, :get_old_state, :get_user_data, :get_user_cache)
    eval(quote
        @inline $op(cb::CoupledCellBuffer) = $op(cb.primary)
    end)
end

get_coupled_buffers(cb::CoupledCellBuffer) = cb.partner_buffers

set_time_increment!(cb::CoupledCellBuffer, Δt) = set_time_increment!(cb.primary, Δt)

for op = (:celldofs, :getcoordinates, :getfieldnames, :cellid)
    eval(quote
        Ferrite.$op(cb::CoupledCellBuffer, args...) = Ferrite.$op(cb.primary, args...)
    end)
end
Ferrite.dof_range(cb::CoupledCellBuffer, name::Symbol) = Ferrite.dof_range(cb.primary, name)

# `reinit_buffer!(cb::CoupledCellBuffer, sim::CoupledSimulation, cellnum::Int)` is defined in
# Coupling.jl, after `CoupledSimulation` exists (this file is included before Coupling.jl).

function reinit_partners!(buffers::NamedTuple, sims::NamedTuple, cellnum::Int)
    map((b, s) -> (reinit_buffer!(b, s, cellnum); nothing), buffers, sims)
    return nothing
end

function _replace_material_with(cb::CoupledCellBuffer, new_material)
    new_primary = _replace_material_with(cb.primary, new_material)
    return CoupledCellBuffer(new_primary, cb.partner_buffers)
end
