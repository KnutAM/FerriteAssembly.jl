using Ferrite.CollectionsOfViews: ArrayOfVectorViews

"""
    QuadratureEvaluator{VT, [type]}(db::Union{DomainBuffer, DomainBuffers})

Overload the function [`eval_quadpoints_cell!`](@ref) to evaluate an expression returning 
the type `VT` for each quadrature point in the domain(s). This is particularly useful in
combination with `Ferrite.L2Projector`. Note that currently only cell domains are supported.
"""
struct QuadratureEvaluator{VT, QEType}
    data::ArrayOfVectorViews{VT}
    qe_type::QEType
    function QuadratureEvaluator(data::ArrayOfVectorViews{VT}, qe_type::Union{Val, Function}) where VT 
        return new{VT, typeof(qe_type)}(data, qe_type)
    end
end
QuadratureEvaluator(data::ArrayOfVectorViews, qe_type::Symbol) = QuadratureEvaluator(data, Val(qe_type))

function QuadratureEvaluator{VT}(domainbuffer::AbstractDomainBuffer, qe_type::Symbol) where {VT}
    cv = get_values(get_itembuffer(domainbuffer))
    @assert cv isa CellValues # For now, only CellValues supported
    nqp = getnquadpoints(cv)
    ncells = getncells(get_dofhandler(domainbuffer).grid)
    data = Vector{VT}(undef, nqp * ncells)
    indices = [1 + i * nqp for i in 0:ncells]
    return QuadratureEvaluator(ArrayOfVectorViews(indices, data, LinearIndices((ncells,))), qe_type)
end

function QuadratureEvaluator{VT}(domainbuffers::DomainBuffers, qe_type::Symbol) where {VT}
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
    return QuadratureEvaluator(ArrayOfVectorViews(indices, data, LinearIndices((ncells,))), qe_type)
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

function eval_quadpoints_cell!(vals::AbstractVector, f::Function, cell_state::AbstractVector, ae, material, cv::CellValues, cellbuffer)
    for q_point in 1:getnquadpoints(cv)
        u  = function_value(v, q_point, ae)
        ∇u = function_gradient(v, q_point, ae)
        vals[q_point] = f(material, u, ∇u, cell_state[q_point])
    end
end

function eval_quadpoints_cell!(vals::AbstractVector, f::Function, cell_state::AbstractVector, ae, material, cv::NamedTuple, cellbuffer)
    for q_point in 1:getnquadpoints(cv)
        u = NamedTuple{keys(fv)}(map((k,v) -> function_value(v, q_point, ae, dof_range(facetbuffer, k)), keys(fv), values(fv)))
        ∇u = NamedTuple{keys(fv)}(map((k,v) -> function_gradient(v, q_point, ae, dof_range(facetbuffer, k)), keys(fv), values(fv)))
        vals[q_point] = f(material, u, ∇u, cell_state[q_point])
    end
end
