mutable struct ElementResidual{S,M,CV,B<:CellBuffer} <: Function
    state::S
    material::M
    cellvalues::CV
    buffer::B
end

function update_element_residual!(er::ElementResidual, state, material, cellvalues, buffer)
    er.state = state
    er.material = material
    er.cellvalues = cellvalues
    er.buffer = buffer 
end

(er::ElementResidual)(re,ae) = element_residual!(re, er.state, ae, er.material, er.cellvalues, er.buffer)

function create_jacobian_config(er::ElementResidual)
    re = get_re(er.buffer)
    ae = get_ae(er.buffer)
    # Setting Chunk explicitly to solve https://github.com/KnutAM/FerriteAssembly.jl/issues/9
    # TODO: This is no longer necessary, but should be benchmarked or at least configurable how
    # to do this. 
    return ForwardDiff.JacobianConfig(er, re, ae, ForwardDiff.Chunk{length(ae)}())
end

struct AutoDiffCellBuffer{CB<:CellBuffer,ER<:ElementResidual,JC} <: AbstractCellBuffer
    cb::CB
    er::ER
    cfg::JC # JacobianConfig
end

include("autodiff_unwrap.jl") # Experimental feature, include to remove large docstring from src here. 

"""
    AutoDiffCellBuffer(cb::CellBuffer)

"""
function AutoDiffCellBuffer(cb::CellBuffer)
    cellstate = deepcopy(get_old_state(cb)) # to be safe, copy shouldn't be required. 
    material = unwrap_material_for_ad(get_material(cb))
    cellvalues = get_values(cb)
    er = ElementResidual(cellstate, material, cellvalues, cb)
    cfg = create_jacobian_config(er)
    return AutoDiffCellBuffer(cb, er, cfg)
end

for op = (:get_Ke, :get_re, :get_ae, :get_material, :get_values, :get_time_increment, 
        :get_aeold, :get_state, :get_old_state, :get_user_data, :get_user_cache, :get_coupled_buffers)
    eval(quote
        @inline $op(cb::AutoDiffCellBuffer) = $op(cb.cb)
    end)
end
reinit_buffer!(cb::AutoDiffCellBuffer, args...; kwargs...) = reinit_buffer!(cb.cb, args...; kwargs...)

set_time_increment!(c::AutoDiffCellBuffer, Δt) = set_time_increment!(c.cb, Δt)

function _replace_material_with(ad_cb::AutoDiffCellBuffer{CB}, new_material) where CB
    cb = setproperties(ad_cb.cb; material = new_material)
    if isa(cb, CB) # If type didn't change, no need to recalculate autodiff buffers
        return setproperties(ad_cb; cb)
    else
        return AutoDiffCellBuffer(cb)
    end
end

function couple_buffers(cb::AutoDiffCellBuffer; kwargs...)
    return AutoDiffCellBuffer(couple_buffers(cb.cb; kwargs...))
end

function create_local(c::AutoDiffCellBuffer)
    cb = create_local(c.cb)
    AutoDiffCellBuffer(cb, deepcopy(c.er), deepcopy(c.cfg))
end

for op = (:celldofs, :getcoordinates, :getfieldnames, :cellid)
    eval(quote
        Ferrite.$op(cb::AutoDiffCellBuffer, args...) = Ferrite.$op(cb.cb, args...)
    end)
end
Ferrite.dof_range(cb::AutoDiffCellBuffer, name::Symbol) = Ferrite.dof_range(cb.cb, name)

function element_routine_ad!(Ke, re, state, ae, material, cellvalues, ad_buffer::AutoDiffCellBuffer)
    buffer = ad_buffer.cb
    element_residual = ad_buffer.er
    cfg = ad_buffer.cfg
    update_element_residual!(element_residual, state, material, cellvalues, buffer)
    try
        ForwardDiff.jacobian!(Ke, element_residual, re, ae, cfg)
    catch e
        ad_error(e, Ke, re, state, ae, material, cellvalues, ad_buffer)
    end
end

# Standard method if no AutoDiffCellBuffer is defined. Should be no need to use, but good to keep for 
# benchmarks if desired. 
function element_routine_ad!(Ke, re, state, ae, material, cellvalues, buffer::CellBuffer)
    rf!(re_, ae_) = element_residual!(re_, state, ae_, material, cellvalues, buffer)
    try
        # Setting Chunk explicitly to solve https://github.com/KnutAM/FerriteAssembly.jl/issues/9
        # TODO: No longer required...
        cfg = ForwardDiff.JacobianConfig(rf!, re, ae, ForwardDiff.Chunk{length(ae)}())
        ForwardDiff.jacobian!(Ke, rf!, re, ae, cfg)
    catch e
        ad_error(e, Ke, re, state, ae, material, cellvalues, buffer)
    end
end

function ad_error(e, Ke, re, state, ae, material, cellvalues, buffer)
    if isa(e, MethodError)
        if e.f === element_residual!
            println(stderr, "Could not find a correctly defined `element_residual!` method")
            showerror(stderr, MethodError(element_routine!, 
                (Ke, re, state, ae, material, cellvalues, buffer))
                ); println()
            println(stderr, "Tried to do automatic differentiation to get Ke,")
            println(stderr, "but could not find a correctly defined `element_residual!` either")
            showerror(stderr, e); println()
            throw(ErrorException("Did not find correctly defined element_routine! or element_residual!"))
        else
            println(stderr, "If the following error is related to converting objects with `ForwardDiff.Dual`s")
            println(stderr, "entries into objects with regular numbers, please consult the docs of `element_residual!`")
        end
    end
    throw(e)
end