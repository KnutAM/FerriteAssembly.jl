
function create_local(a::Ferrite.AbstractAssembler)
    start_assemble(Ferrite.get_matrix(a), Ferrite.get_vector(a); fillzero=false)
end

scatter!(::Ferrite.AbstractAssembler, ::Ferrite.AbstractAssembler) = nothing
gather!( ::Ferrite.AbstractAssembler, ::Ferrite.AbstractAssembler) = nothing

can_thread(::Ferrite.AbstractAssembler) = true

function work_single_cell!(assembler::Ferrite.AbstractAssembler, cellbuffer)
    Ke = get_Ke(cellbuffer)
    re = get_re(cellbuffer)
    ae = get_ae(cellbuffer)
    cell_state = get_state(cellbuffer)
    material = get_material(cellbuffer)
    cellvalues = get_values(cellbuffer)
    element_routine!(Ke, re, cell_state, ae, material, cellvalues, cellbuffer)
    assemble!(assembler, celldofs(cellbuffer), Ke, re)
end

function work_single_facet!(assembler::Ferrite.AbstractAssembler, facetbuffer)
    Ke = get_Ke(facetbuffer)
    re = get_re(facetbuffer)
    ae = get_ae(facetbuffer)
    material = get_material(facetbuffer)
    cellvalues = get_values(facetbuffer)
    facet_routine!(Ke, re, ae, material, cellvalues, facetbuffer)
    assemble!(assembler, celldofs(facetbuffer), Ke, re)
end

# Assemble only the residual
"""
    ReAssembler(r; fillzero=true, scaling=NoScaling())
    
`ReAssembler` will call `element_residual!` and assemble `re` into the system 
vector `r`. `fillzero=true` implies that `r` is zeroed upon construction of 
`ReAssembler`. Furthermore, remaining keyword arguments can enable the following 
features:
- `scaling`: Calculate a scaling measure locally at the residual level, see e.g., 
  [`ElementResidualScaling`](@ref). 
  
[`reset_scaling!`](@ref) is called when creating `ReAssembler`
"""
struct ReAssembler{R<:AbstractVector, SC}
    r::R
    scaling::SC
end
function ReAssembler(r::AbstractVector{T}; scaling=NoScaling(), fillzero=true) where T
    fillzero && fill!(r, zero(T))
    reset_scaling!(scaling)
    return ReAssembler(r, scaling)
end

# TaskLocals interface:
function create_local(ra::ReAssembler)
    scaling = create_local(ra.scaling)
    return ReAssembler(ra.r, scaling)
end
scatter!(task::ReAssembler, base::ReAssembler) = scatter!(task.scaling, base.scaling)
gather!(base::ReAssembler, task::ReAssembler) = gather!(base.scaling, task.scaling)

can_thread(::ReAssembler) = true

function work_single_cell!(assembler::ReAssembler, cellbuffer)
    re = get_re(cellbuffer)
    ae = get_ae(cellbuffer)
    cell_state = get_state(cellbuffer)
    material = get_material(cellbuffer)
    cellvalues = get_values(cellbuffer)
    element_residual!(re, cell_state, ae, material, cellvalues, cellbuffer)
    update_scaling!(assembler.scaling, re, cellbuffer)
    assemble!(assembler.r, celldofs(cellbuffer), re)
end

function work_single_facet!(assembler::ReAssembler, facetbuffer)
    re = get_re(facetbuffer)
    ae = get_ae(facetbuffer)
    material = get_material(facetbuffer)
    facetvalues = get_values(facetbuffer)
    facet_residual!(re, ae, material, facetvalues, facetbuffer)
    update_scaling!(assembler.scaling, re, facetbuffer)
    assemble!(assembler.r, celldofs(facetbuffer), re)
end

"""
    KeReAssembler(K, r; fillzero=true, kwargs...)
    KeReAssembler(a::Ferrite.AbstractAssembler; kwargs...)

The default `KeReAssembler` works just like a `Ferrite.AbstractAssembler`:
It will call `element_routine!` and assemble `Ke` and `re` into the global system 
matrix `K` and vector `r`. However, it comes with the additional possible features,
that are controllable via the keyword arguments:
- `ch::Union{Nothing,ConstraintHandler}=nothing`: If a `ConstraintHandler` is given,
  local applications of constraints will be applied, using `Ferrite.apply_assemble`.
- `apply_zero`: Required if a constraint handler is given, and forwarded to `Ferrite.apply_assemble`. 
- `scaling`: Calculate a scaling measure locally at the residual level, see e.g., 
  [`ElementResidualScaling`](@ref)

`a=Ferrite.start_assemble(K, r; fillzero=fillzero)` is passed to the second definition 
if a matrix, `K`, and vector, `r`, are given as input.

[`reset_scaling!`](@ref) is called when creating `KeReAssembler`
"""
struct KeReAssembler{A<:Ferrite.AbstractAssembler, CH, SC}
    a::A
    ch::CH # Union{Nothing, Ferrite.ConstraintHandler}
    apply_zero::Bool
    scaling::SC
end
function KeReAssembler(K::AbstractMatrix, r::AbstractVector; fillzero=true, kwargs...)
    a = start_assemble(K, r; fillzero=fillzero)
    return KeReAssembler(a; kwargs...)
end
function KeReAssembler(a::Ferrite.AbstractAssembler; apply_zero=nothing, ch=nothing, scaling=NoScaling())
    reset_scaling!(scaling)
    if !isnothing(ch) && isnothing(apply_zero)
         throw(ArgumentError("apply_zero must be specified when `ch` is given"))
    end
    _apply_zero = isnothing(apply_zero) ? false : apply_zero
    KeReAssembler(a, ch, _apply_zero, scaling)
end

# TaskLocals interface:
function create_local(ra::KeReAssembler)
    a = create_local(ra.a)
    scaling = create_local(ra.scaling)
    return KeReAssembler(a, ra.ch, ra.apply_zero, scaling)
end
function scatter!(task::KeReAssembler, base::KeReAssembler)
    scatter!(task.a, base.a)
    scatter!(task.scaling, base.scaling)
end
function gather!(base::KeReAssembler, task::KeReAssembler)
    gather!(base.a, task.a)
    gather!(base.scaling, task.scaling)
end

can_thread(::KeReAssembler) = true

# assemble! routines
# # No constraint handler - no local application of constraints
function assemble_contributions!(assembler::KeReAssembler{<:Any,Nothing}, buffer::AbstractItemBuffer)
    Ke = get_Ke(buffer)
    re = get_re(buffer)
    update_scaling!(assembler.scaling, re, buffer)
    assemble!(assembler.a, celldofs(buffer), Ke, re)
end
# Constraint handler - application of constraints
function assemble_contributions!(assembler::KeReAssembler{<:Any,<:ConstraintHandler}, buffer::AbstractItemBuffer)
    Ke = get_Ke(buffer)
    re = get_re(buffer)
    update_scaling!(assembler.scaling, re, buffer)
    apply_assemble!(assembler.a, assembler.ch, celldofs(buffer), Ke, re; apply_zero=assembler.apply_zero)
end

function work_single_cell!(assembler::KeReAssembler, cellbuffer)
    Ke = get_Ke(cellbuffer)
    re = get_re(cellbuffer)
    ae = get_ae(cellbuffer)
    cell_state = get_state(cellbuffer)
    material = get_material(cellbuffer)
    cellvalues = get_values(cellbuffer)
    element_routine!(Ke, re, cell_state, ae, material, cellvalues, cellbuffer)
    assemble_contributions!(assembler, cellbuffer)
end

function work_single_facet!(assembler::KeReAssembler, facetbuffer)
    Ke = get_Ke(facetbuffer)
    re = get_re(facetbuffer)
    ae = get_ae(facetbuffer)
    material = get_material(facetbuffer)
    cellvalues = get_values(facetbuffer)
    facet_routine!(Ke, re, ae, material, cellvalues, facetbuffer)
    assemble_contributions!(assembler, facetbuffer)
end
