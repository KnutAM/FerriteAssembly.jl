using Ferrite.CollectionsOfViews: ArrayOfVectorViews

"""
    QuadPointEvaluator{VT}(db::Union{DomainBuffer, DomainBuffers}, f::Function)

Create a `QuadPointEvaluator` for the domain(s) `db`, for which the function
`f(material, u, ∇u, qp_state)` will be used to calculate a value of type `VT` in 
each quadrature point. A `QuadPointEvaluator` is often used in combination with 
`Ferrite`'s `L2Projector`. Note that currently only cell domains are supported.

Currently, the stored data should be accessed directly via the `data` field. 
`data` is indexed by the cell number, giving an `AbstractVector` of the values 
for each quadrature point in that cell. A more proper API may be introduced in the future. 

    QuadPointEvaluator{VT}(db::Union{DomainBuffer, DomainBuffers}, qe_type::Symbol)

A more flexible version of a `QuadPointEvaluator` where it is required to overload
[`eval_quadpoints_cell!`](@ref) to calculate the values to be stored for each quadrature point. 
"""
struct QuadPointEvaluator{VT, QEType}
    data::ArrayOfVectorViews{VT}
    qe_type::QEType
    function QuadPointEvaluator(data::ArrayOfVectorViews{VT}, qe_type::Union{Val, Function}) where VT 
        return new{VT, typeof(qe_type)}(data, qe_type)
    end
end
QuadPointEvaluator(data::ArrayOfVectorViews, qe_type::Symbol) = QuadPointEvaluator(data, Val(qe_type))

function QuadPointEvaluator{VT}(domainbuffer::AbstractDomainBuffer, qe_type::Union{Symbol, Function}) where {VT}
    cv = get_values(get_itembuffer(domainbuffer))
    nqp = if cv isa Ferrite.AbstractCellValues
        getnquadpoints(cv)
    elseif cv isa NamedTuple
        tmp = getnquadpoints.(values(cv))
        @assert allequal(tmp)
        tmp[1]
    else
        error("Only CellValues are supported")
    end
    ncells = getncells(get_dofhandler(domainbuffer).grid)
    data = Vector{VT}(undef, nqp * ncells)
    indices = [1 + i * nqp for i in 0:ncells]
    return QuadPointEvaluator(ArrayOfVectorViews(indices, data, LinearIndices((ncells,))), qe_type)
end

function QuadPointEvaluator{VT}(domainbuffers::DomainBuffers, qe_type::Union{Symbol, Function}) where {VT}
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
    return QuadPointEvaluator(ArrayOfVectorViews(indices, data, LinearIndices((ncells,))), qe_type)
end

can_thread(::QuadPointEvaluator) = true 
skip_this_domain(::QuadPointEvaluator, ::String) = false

"""
    eval_quadpoints_cell!(qp_values::AbstractVector, ::Val{type}, cell_state, ae, material, cv, cellbuffer)

Function to be overloaded for `Val{type}`, potentially in combination with `material`.
Will be called when using `QuadPointEvaluator` to add evaluated values for each quadrature point in 
`cv` to qp_values. If no `type` is given to `QuadPointEvaluator`, `type=:default` will be used. 
"""
function eval_quadpoints_cell! end

function work_single_cell!(qe::QuadPointEvaluator, cellbuffer)
    cv = get_values(cellbuffer)
    m = get_material(cellbuffer)
    cell_state = get_state(cellbuffer)
    ae = get_ae(cellbuffer)
    eval_quadpoints_cell!(qe.data[cellid(cellbuffer)], qe.qe_type, cell_state, ae, m, cv, cellbuffer)
end

function eval_quadpoints_cell!(vals::AbstractVector, f::Function, cell_state::AbstractVector, ae, material, cv::AbstractCellValues, cellbuffer)
    for q_point in 1:getnquadpoints(cv)
        u  = function_value(cv, q_point, ae)
        ∇u = function_gradient(cv, q_point, ae)
        vals[q_point] = f(material, u, ∇u, cell_state[q_point])
    end
end

function eval_quadpoints_cell!(vals::AbstractVector, f::Function, cell_state::AbstractVector, ae, material, cv::NamedTuple, cellbuffer)
    for q_point in 1:getnquadpoints(first(cv))
        u = NamedTuple{keys(cv)}(map((k,v) -> function_value(v, q_point, ae, dof_range(cellbuffer, k)), keys(cv), values(cv)))
        ∇u = NamedTuple{keys(cv)}(map((k,v) -> function_gradient(v, q_point, ae, dof_range(cellbuffer, k)), keys(cv), values(cv)))
        vals[q_point] = f(material, u, ∇u, cell_state[q_point])
    end
end
