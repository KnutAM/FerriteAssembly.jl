mutable struct ElementResidual{S,M,CV,DH_FH,T,B<:CellBuffer} <: Function
    state::S
    material::M
    cellvalues::CV
    dh_fh::DH_FH
    Δt::T 
    buffer::B
end

function update_element_residual!(er::ElementResidual, state, material, cellvalues, dh_fh, Δt, buffer)
    er.state = state
    er.material = material
    er.cellvalues = cellvalues
    er.dh_fh = dh_fh
    er.Δt = isa(Δt,Nothing) ? zero(er.Δt) : Δt
    er.buffer = buffer 
end

(er::ElementResidual)(re,ae) = element_residual!(re, er.state, ae, er.material, er.cellvalues, er.dh_fh, er.Δt, er.buffer)

function create_jacobian_config(er::ElementResidual)
    re = get_re(er.buffer)
    ae = get_ae(er.buffer)
    # Setting Chunk explicitly to solve https://github.com/KnutAM/FerriteAssembly.jl/issues/9
    return ForwardDiff.JacobianConfig(er, re, ae, ForwardDiff.Chunk{length(ae)}())
end

struct AutoDiffCellBuffer{CB<:CellBuffer,ER<:ElementResidual,JC} <: AbstractCellBuffer
    cb::CB
    er::ER
    cfg::JC # JacobianConfig
end

# Mini_helpers
extract_first_cellstate(s) = first(s)
extract_first_cellstate(s::Dict) = s[first(keys(s))]    # first(s) returns Pair(key,val)

"""
    AutoDiffCellBuffer(cb::CellBuffer, cellstates, dh_fh::Union{DofHandler,FieldHandler})

`cellstates` should be the states for the relevant cells in `dh_fh`. 
An `AutoDiffCellBuffer` wraps a `CellBuffer` and behaves the same way 
wrt defined functions, i.e. all `get*` functions defined for `CellBuffer`
are also defined for `AutoDiffCellBuffer`. 
When used for automatic differentiation, only the wrapped `CellBuffer` 
is passed to `element_residual!`. 
As a backup solution if direct field access is used, the function 
[`getCellBuffer`](@ref) returns the `CellBuffer` for both a `CellBuffer` 
and an `AutoDiffCellBuffer`. 

Note that this constructor is normally not used, 
but called from the [`setup_ad_cellbuffer`](@ref) function.
"""
function AutoDiffCellBuffer(cb::CellBuffer, cellstates, dh_fh::Union{DofHandler,FieldHandler})
    cellstate = extract_first_cellstate(cellstates)
    Δt=zero(Float64)
    material = get_material(cb)
    cellvalues = get_cellvalues(cb)
    er = ElementResidual(cellstate, material, cellvalues, dh_fh, Δt, cb)
    cfg = create_jacobian_config(er)
    return AutoDiffCellBuffer(cb, er, cfg)
end

"""
    setup_ad_cellbuffer(states, dh::AbstractDofHandler, args...)

Create a special cell buffer that improves performance when using 
automatic differentiation to calculate the element stiffness. The exact same 
inputs as for [`setup_cellbuffer`](@ref) should be used, except that the 
`states` variable that is passed to `doassemble!` should be prepended to the 
argument list.  
"""
function setup_ad_cellbuffer end

# General 
setup_ad_cellbuffer(cb::CellBuffer, args...) = AutoDiffCellBuffer(cb, args...)

# DofHandler
function setup_ad_cellbuffer(states, dh::DofHandler, args...)
    AutoDiffCellBuffer(setup_cellbuffer(dh, args...), states, dh)
end
# MixedDofHandler
function setup_ad_cellbuffer(states::Tuple, dh::MixedDofHandler, args...)
    return setup_ad_cellbuffer(setup_cellbuffer(dh, args...), states, dh)
end
function setup_ad_cellbuffer(cbs::Tuple, states::Tuple, dh::MixedDofHandler)
    numfh = length(dh.fieldhandlers)
    numfh==length(states) || error("length(states)=$(length(states))!=$numfh=length(fieldhandlers)")
    return ntuple(i->AutoDiffCellBuffer(cbs[i], states[i], dh.fieldhandlers[i]), numfh)
end
# [Mixed]DofHandler + mixed materials
function setup_ad_cellbuffer(states::Dict{String}, dh::DofHandler, args...)
    return _setup_dict_ad_cellbuffer(states, dh, args...)
end
function setup_ad_cellbuffer(states::Dict{String}, dh::MixedDofHandler, args...)
    return _setup_dict_ad_cellbuffer(states, dh, args...)
end
# Need special case to avoid type ambiguity
function _setup_dict_ad_cellbuffer(states::Dict{String}, dh, cv, mtrls::Dict{String}, args...)
    keys(states) == keys(mtrls) || error("states::Dict and mtrls::Dict must have the same keys")
    cbs = setup_cellbuffer(dh, cv, mtrls, args...)
    return Dict(key=>setup_ad_cellbuffer(cbs[key], states[key], dh) for key in keys(states))
end

"""
    getCellBuffer(b::Union{CellBuffer,AutoDiffCellBuffer})

Return the `CellBuffer` if `isa(b,CellBuffer)`, else, return the `CellBuffer`
wrapped by `b::AutoDiffCellBuffer`, i.e. `b.cb`
"""
@inline getCellBuffer(b::AutoDiffCellBuffer) = b.cb
@inline getCellBuffer(b::CellBuffer) = b

for op = (:get_Ke, :get_re, :get_ae, :get_material, :get_cellvalues, :get_aeold, :get_load, :get_cache)
    eval(quote
        $op(cb::AutoDiffCellBuffer) = $op(cb.cb)
    end)
end
for op = (:celldofs, :getcoordinates, :reinit!)
    eval(quote
        Ferrite.$op(cb::AutoDiffCellBuffer, args...) = Ferrite.$op(cb.cb, args...)
    end)
end

function element_routine_ad!(Ke, re, state, ae, material, cellvalues, dh_fh, Δt, ad_buffer::AutoDiffCellBuffer)
    buffer = getCellBuffer(ad_buffer)
    element_residual = ad_buffer.er
    cfg = ad_buffer.cfg
    update_element_residual!(element_residual, state, material, cellvalues, dh_fh, Δt, buffer)
    try
        ForwardDiff.jacobian!(Ke, element_residual, re, ae, cfg)
    catch e
        ad_error(e, Ke, re, state, ae, material, cellvalues, dh_fh, Δt, ad_buffer)
    end
end

# Standard method if no AutoDiffCellBuffer is defined. Should be no need to use, but good to keep for 
# benchmarks if desired. 
function element_routine_ad!(Ke, re, state, ae, material, cellvalues, dh_fh, Δt, buffer::CellBuffer)
    rf!(re_, ae_) = element_residual!(re_, state, ae_, material, cellvalues, dh_fh, Δt, buffer)
    try
        # Setting Chunk explicitly to solve https://github.com/KnutAM/FerriteAssembly.jl/issues/9
        cfg = ForwardDiff.JacobianConfig(rf!, re, ae, ForwardDiff.Chunk{length(ae)}())
        ForwardDiff.jacobian!(Ke, rf!, re, ae, cfg)
    catch e
        ad_error(e, Ke, re, state, ae, material, cellvalues, dh_fh, Δt, buffer)
    end
end

function ad_error(e, Ke, re, state, ae, material, cellvalues, dh_fh, Δt, buffer)
    if isa(e, MethodError)
        if e.f === element_residual!
            println("Could not find a correctly defined `element_residual!` method")
            showerror(Base.stderr, MethodError(element_routine!, 
                (Ke, re, state, ae, material, cellvalues, dh_fh, Δt, buffer))
                ); println()
            println("Tried to do automatic differentiation to get Ke,")
            println("but could not find a correctly defined `element_residual!` either")
            showerror(Base.stderr, e); println()
            throw(ErrorException("Did not find correctly defined element_routine! or element_residual!"))
        else
            println("If the following error is related to converting objects with `ForwardDiff.Dual`s")
            println("entries into objects with regular numbers, please consult the docs of `element_residual!`")
        end
    end
    rethrow()
end