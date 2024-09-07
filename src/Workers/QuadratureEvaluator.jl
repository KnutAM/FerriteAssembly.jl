using Ferrite.CollectionsOfViews: ArrayOfVectorViews

"""
    QuadratureEvaluator{VT, [type]}(db::Union{DomainBuffer, DomainBuffers})

Overload the function [`eval_quadpoints_cell!`](@ref) to evaluate an expression returning 
the type `VT` for each quadrature point in the domain(s). This is particularly useful in
combination with `Ferrite.L2Projector`. Note that currently only cell domains are supported.
"""
struct QuadratureEvaluator{VT, type}
    data::ArrayOfVectorViews{VT}
end
QuadratureEvaluator(data::ArrayOfVectorViews{VT}) where VT = QuadratureEvaluator{VT, :default}(data)

function QuadratureEvaluator{VT}(domainbuffer::AbstractDomainBuffer) where {VT}
    cv = get_values(get_itembuffer(domainbuffer))
    @assert cv isa CellValues # For now, only CellValues supported
    nqp = getnquadpoints(cv)
    ncells = getncells(get_dofhandler(domainbuffer).grid)
    data = Vector{VT}(undef, nqp * ncells)
    indices = [1 + i * nqp for i in 0:ncells]
    return QuadratureEvaluator(ArrayOfVectorViews(indices, data, LinearIndices((ncells,))))
end

function QuadratureEvaluator{VT}(domainbuffers::DomainBuffers) where {VT}
    ncells = getncells(get_dofhandler(first(values(domainbuffers))).grid)
    nqps = Vector{Int}(undef, ncells)
    nqp_total = 0
    for (_, db) in domainbuffers
        cv = get_values(get_itembuffer(db))
        @assert cv isa CellValues # For now, only CellValues supported
        nqp = getnquadpoints(cv)
        set = getset(db)
        map(i -> (nqps[i] = nqp), set)
        nqp_total += nqp * length(set)
    end

    data = Vector{VT}(undef, nqp_total)
    indices = Vector{Int}(undef, ncells + 1)
    indices[1] = 1
    for (i, nqp) in enumerate(nqps)
        indices[i + 1] = indices[i] + nqp
    end
    return QuadratureEvaluator(ArrayOfVectorViews(indices, data, LinearIndices((ncells,))))
end

can_thread(::QuadratureEvaluator) = true 
skip_this_domain(::QuadratureEvaluator, ::String) = false

"""
    eval_quadpoints_cell!(qp_values::AbstractVector, ::Val{type}, cell_state, ae, material, cv, cellbuffer)

Function to be overloaded for `Val{type}`, potentially in combination with `material`.
Will be called when using `QuadratureEvaluator` to add evaluated values for each quadrature point in 
`cv` to qp_values. If no `type` is given to `QuadratureEvaluator`, `type=:default` will be used. 
"""
function eval_quadpoints_cell! end

function work_single_cell!(qe::QuadratureEvaluator{<:Any, type}, cellbuffer) where {type}
    cv = get_values(cellbuffer)
    m = get_material(cellbuffer)
    cell_state = get_state(cellbuffer)
    ae = get_ae(cellbuffer)
    eval_quadpoints_cell!(qe.data[getcellid(cellbuffer)], Val(type), cell_state, ae, m, cv, cellbuffer)
end
