"""
    CoupledCellBuffer(primary::CellBuffer, partner_buffers::NamedTuple, partner_sims::NamedTuple)

Wraps a reader's own `primary::CellBuffer` together with references to its coupling
partners' plain `CellBuffer`s (`partner_buffers`, returned by [`get_coupled_buffer`](@ref))
and the partner `Simulation`s required to reinitialize them (`partner_sims`, internal only).

`partner_buffers` are never themselves `CoupledCellBuffer`s or `AutoDiffCellBuffer`s:
coupling is direct and nonrecursive, so a partner's own coupling (if any) is not exposed.

Constructed once per task at [`CoupledSimulations`](@ref) setup time; never rebuilt per cell
or per `work!` call.
"""
struct CoupledCellBuffer{CB<:CellBuffer, PB<:NamedTuple, PS<:NamedTuple} <: AbstractCellBuffer
    primary::CB
    partner_buffers::PB
    partner_sims::PS
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

"""
    reinit_buffer!(cb::CoupledCellBuffer, sim::Simulation, cellnum::Int)

Reinitialize the reader's own `cb.primary` against `sim`, then reinitialize each
partner buffer in `cb.partner_buffers` against its own partner `Simulation` stored in
`cb.partner_sims`. Partner reinitialization does not recurse: partner buffers are plain
`CellBuffer`s, so no further coupling initialization happens.
"""
function reinit_buffer!(cb::CoupledCellBuffer, sim::Simulation, cellnum::Int)
    reinit_buffer!(cb.primary, sim, cellnum)
    reinit_partners!(cb.partner_buffers, cb.partner_sims, cellnum)
    return nothing
end

function reinit_partners!(buffers::NamedTuple, sims::NamedTuple, cellnum::Int)
    map((b, s) -> (reinit_buffer!(b, s, cellnum); nothing), buffers, sims)
    return nothing
end

function _replace_material_with(cb::CoupledCellBuffer, new_material)
    new_primary = _replace_material_with(cb.primary, new_material)
    return CoupledCellBuffer(new_primary, cb.partner_buffers, cb.partner_sims)
end
