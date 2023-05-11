# Assemblers need to the following:
# `TaskLocals` interface
# `assemble_cell_reinited!(assembler, state, cellbuffer)`

# Ferrite.AbstractSparseAssembler
const FerriteSparseAssemblers = Union{Ferrite.AssemblerSparsityPattern, Ferrite.AssemblerSymmetricSparsityPattern}
create_local(a::FerriteSparseAssemblers) = start_assemble(a.K, a.f; fillzero=false)
scatter!(::FerriteSparseAssemblers, ::FerriteSparseAssemblers) = nothing
gather!( ::FerriteSparseAssemblers, ::FerriteSparseAssemblers) = nothing

"""
    assemble_cell_reinited!(assembler::Ferrite.AbstractSparseAssembler, cell_state, cellbuffer)

Calculate the local element stiffness, `Ke`, and the local residual, `re`, with
`element_routine!` for the cell described by the reinitialized `cellbuffer`.
Assemble these into the global matrix and vector in `assembler`. 
"""
function assemble_cell_reinited!(assembler::Ferrite.AbstractSparseAssembler, cell_state, cellbuffer)
    Ke = get_Ke(cellbuffer)
    re = get_re(cellbuffer)
    ae = get_ae(cellbuffer)
    material = get_material(cellbuffer)
    cellvalues = get_cellvalues(cellbuffer)
    element_routine!(Ke, re, cell_state, ae, material, cellvalues, cellbuffer)
    assemble!(assembler, celldofs(cellbuffer), Ke, re)
end

# Assemble only the residual
"""
    ReAssembler(r; fillzero=true, scaling=NoScaling())
    
`ReAssembler` will call `element_residual!` and assemble `re` into the system 
vector `r`. `fillzero=true` implies that `r` is zeroed upon construction of 
`ReAssembler`. Furthermore, remaining keyword arguments can enable the following 
features:
- `scaling`: Calculate a scaling measure locally at the residual level, see e.g., 
  [`ElementResidualScaling`](@ref)
"""
struct ReAssembler{R<:AbstractVector, SC}
    r::R
    scaling::SC
end
function ReAssembler(r::AbstractVector{T}; scaling=NoScaling(), fillzero=true) where T
    fillzero && fill!(r, zero(T))
    return ReAssembler(r, scaling)
end

# TaskLocals interface:
function create_local(ra::ReAssembler)
    scaling = create_local(ra.scaling)
    return ReAssembler(ra.r, scaling)
end
scatter!(task::ReAssembler, base::ReAssembler) = scatter!(task.scaling, base.scaling)
gather!(base::ReAssembler, task::ReAssembler) = gather!(base.scaling, task.scaling)

# assemble! routine
function assemble_contributions!(ra::ReAssembler, buffer::AbstractCellBuffer)
    re = get_re(buffer)
    update_scaling!(ra.scaling, re, buffer)
    assemble!(ra.r, celldofs(buffer), re)
end

"""
    assemble_cell_reinited!(assembler::ReAssembler, cell_state, cellbuffer)

Calculate the local residual, `re`, with `element_residual!` for the cell described 
by the reinitialized `cellbuffer`. Assemble `re` into the global vector in `assembler`. 
In addition, `ReAssembler` supports the extra option
- Residual scaling factor, e.g. `ElementResidualScaling`
"""
function assemble_cell_reinited!(assembler::ReAssembler, cell_state, cellbuffer)
    re = get_re(cellbuffer)
    ae = get_ae(cellbuffer)
    material = get_material(cellbuffer)
    cellvalues = get_cellvalues(cellbuffer)
    element_residual!(re, cell_state, ae, material, cellvalues, cellbuffer)
    
    assemble_contributions!(assembler, cellbuffer)
end

"""
    KeReAssembler(K, r; fillzero=true, kwargs...)
    KeReAssembler(a::Ferrite.AbstractSparseAssembler; kwargs...)

The default `KeReAssembler` works just like a `Ferrite.AbstractSparseAssembler`:
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
"""
struct KeReAssembler{A<:FerriteSparseAssemblers, CH, SC}
    a::A
    ch::CH # Union{Nothing, Ferrite.ConstraintHandler}
    apply_zero::Bool
    scaling::SC
end
function KeReAssembler(K::AbstractMatrix, r::AbstractVector; fillzero=true, kwargs...)
    a = start_assemble(K, r; fillzero=fillzero)
    return KeReAssembler(a; kwargs...)
end
function KeReAssembler(a::FerriteSparseAssemblers; apply_zero=nothing, ch=nothing, scaling=NoScaling())
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

# assemble! routines
# # No constraint handler - no local application of constraints
function assemble_contributions!(assembler::KeReAssembler{<:Any,Nothing}, buffer::AbstractCellBuffer)
    Ke = get_Ke(buffer)
    re = get_re(buffer)
    update_scaling!(assembler.scaling, re, buffer)
    assemble!(assembler.a, celldofs(buffer), Ke, re)
end
# Constraint handler - application of constraints
function assemble_contributions!(assembler::KeReAssembler{<:Any,<:ConstraintHandler}, buffer::AbstractCellBuffer)
    Ke = get_Ke(buffer)
    re = get_re(buffer)
    update_scaling!(assembler.scaling, re, buffer)
    apply_assemble!(assembler.a, assembler.ch, celldofs(buffer), Ke, re; apply_zero=assembler.apply_zero)
end

"""
    assemble_cell_reinited!(assembler::KeReAssembler, cell_state, cellbuffer)

Calculate the local element stiffness, `Ke`, and the local residual, `re`, with
`element_routine!` for the cell described by the reinitialized `cellbuffer`.
Assemble these into the global matrix and vector in `assembler`. 
In addition, `KeReAssembler` supports the extra options
- Local application of constraints: `Ferrite.apply_assemble!`
- Residual scaling factor, e.g. `ElementResidualScaling`
"""
function assemble_cell_reinited!(assembler::KeReAssembler, cell_state, cellbuffer)
    Ke = get_Ke(cellbuffer)
    re = get_re(cellbuffer)
    ae = get_ae(cellbuffer)
    material = get_material(cellbuffer)
    cellvalues = get_cellvalues(cellbuffer)
    element_routine!(Ke, re, cell_state, ae, material, cellvalues, cellbuffer)
    assemble_contributions!(assembler, cellbuffer)
end
