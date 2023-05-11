# Assemblers should support the `TaskLocals` structure

# Ferrite.AbstractSparseAssembler
const FerriteSparseAssemblers = Union{Ferrite.AssemblerSparsityPattern, Ferrite.AssemblerSymmetricSparsityPattern}
create_local(a::FerriteSparseAssemblers) = start_assemble(a.K, a.f; fillzero=false)
scatter!(::FerriteSparseAssemblers, ::FerriteSparseAssemblers) = nothing
gather!( ::FerriteSparseAssemblers, ::FerriteSparseAssemblers) = nothing

# Assemble only the residual
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
function Ferrite.assemble!(ra::ReAssembler, buffer::AbstractCellBuffer)
    re = get_re(buffer)
    update_scaling!(ra.scaling, re, buffer)
    for (i, I) in pairs(celldofs(buffer))
        ra.r[I] += re[i]
    end
end

# Assembler that builds on top of existing assemblers
struct KeReAssembler{A<:FerriteSparseAssemblers, CH, SC}
    a::A
    ch::CH # Union{Nothing, Ferrite.ConstraintHandler}
    apply_zero::Bool
    scaling::SC
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
function Ferrite.assemble!(assembler::KeReAssembler{<:Any,Nothing}, buffer::AbstractCellBuffer)
    Ke = get_Ke(buffer)
    re = get_re(buffer)
    update_scaling!(assembler.scaling, re, buffer)
    assemble!(assembler.a, celldofs(buffer), Ke, re)
end
# Constraint handler - application of constraints
function Ferrite.assemble!(assembler::KeReAssembler{<:Any,<:ConstraintHandler}, buffer::AbstractCellBuffer)
    Ke = get_Ke(buffer)
    re = get_re(buffer)
    update_scaling!(assembler.scaling, re, buffer)
    apply_assemble!(assembler.a, assembler.ch, celldofs(buffer), Ke, re; apply_zero=assembler.apply_zero)
end
