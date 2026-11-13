# Values types with additional base functions, used to test variables attached to fields
# (`AttachAlgebraic` in the `algebraic_variables` of a `DomainSpec`)

# Lift 2d tensors to 3d, padding with zeros
lift3d(v::Number) = v
lift3d(v::Vec{2, T}) where {T} = Vec{3, T}((v[1], v[2], zero(T)))
lift3d(t::Tensor{2, 2, T}) where {T} = Tensor{2, 3, T}((i, j) -> (i ≤ 2 && j ≤ 2) ? t[i, j] : zero(T))
lift3d(t::Tensor{3, 2, T}) where {T} = Tensor{3, 3, T}((i, j, k) -> (i ≤ 2 && j ≤ 2 && k ≤ 2) ? t[i, j, k] : zero(T))

# 2d cell values lifted to 3d, with one additional base function, which has a zero value and the gradient
# e₃ (scalar field) or e₃ ⊗ e₃ (vector field), i.e. a constant out-of-plane gradient (generalized plane strain).
struct OutOfPlaneCellValues{CV <: CellValues} <: Ferrite.AbstractCellValues
    cv::CV
end
Ferrite.getnbasefunctions(v::OutOfPlaneCellValues) = getnbasefunctions(v.cv) + 1
Ferrite.getnquadpoints(v::OutOfPlaneCellValues) = getnquadpoints(v.cv)
Ferrite.getdetJdV(v::OutOfPlaneCellValues, q_point::Int) = getdetJdV(v.cv, q_point)
Ferrite.getngeobasefunctions(v::OutOfPlaneCellValues) = Ferrite.getngeobasefunctions(v.cv)
Ferrite.geometric_value(v::OutOfPlaneCellValues, q_point::Int, i::Int) = Ferrite.geometric_value(v.cv, q_point, i)
function Ferrite.reinit!(v::OutOfPlaneCellValues, cell::Union{Nothing, Ferrite.AbstractCell}, x::AbstractVector{<:Vec})
    return reinit!(v.cv, cell, x)
end
out_of_plane_gradient(::Number) = Vec{3}((0.0, 0.0, 1.0))
out_of_plane_gradient(::Vec) = Tensor{2, 3}((i, j) -> (i == 3 && j == 3) ? 1.0 : 0.0)
function Ferrite.shape_value(v::OutOfPlaneCellValues, q_point::Int, i::Int)
    N = lift3d(shape_value(v.cv, q_point, 1))
    return i ≤ getnbasefunctions(v.cv) ? lift3d(shape_value(v.cv, q_point, i)) : zero(N)
end
function Ferrite.shape_gradient(v::OutOfPlaneCellValues, q_point::Int, i::Int)
    i ≤ getnbasefunctions(v.cv) && return lift3d(shape_gradient(v.cv, q_point, i))
    return out_of_plane_gradient(shape_value(v.cv, q_point, 1))
end
function Ferrite.shape_hessian(v::OutOfPlaneCellValues, q_point::Int, i::Int)
    H = lift3d(shape_hessian(v.cv, q_point, 1))
    return i ≤ getnbasefunctions(v.cv) ? lift3d(shape_hessian(v.cv, q_point, i)) : zero(H)
end

# Vector-valued cell values extended with the modes of the field u = H ⋅ (x - xc) + g, for the
# given components, `H_indices`, of the tensor H and the components, `g_indices`, of the vector g.
# The positions x - xc in the quadrature points are calculated in `reinit!`.
struct AffineCellValues{dim, CV <: CellValues, NH, NG} <: Ferrite.AbstractCellValues
    cv::CV
    xc::Vec{dim, Float64}
    H_indices::NTuple{NH, Tuple{Int, Int}}
    g_indices::NTuple{NG, Int}
    Δx::Vector{Vec{dim, Float64}}
end
function AffineCellValues(cv::CellValues, xc::Vec{dim}, H_indices, g_indices = ()) where {dim}
    return AffineCellValues(cv, Vec{dim, Float64}(xc), H_indices, g_indices, zeros(Vec{dim, Float64}, getnquadpoints(cv)))
end
Ferrite.getnbasefunctions(v::AffineCellValues) = getnbasefunctions(v.cv) + length(v.H_indices) + length(v.g_indices)
Ferrite.getnquadpoints(v::AffineCellValues) = getnquadpoints(v.cv)
Ferrite.getdetJdV(v::AffineCellValues, q_point::Int) = getdetJdV(v.cv, q_point)
Ferrite.getngeobasefunctions(v::AffineCellValues) = Ferrite.getngeobasefunctions(v.cv)
Ferrite.geometric_value(v::AffineCellValues, q_point::Int, i::Int) = Ferrite.geometric_value(v.cv, q_point, i)
function Ferrite.reinit!(v::AffineCellValues, cell::Union{Nothing, Ferrite.AbstractCell}, x::AbstractVector{<:Vec})
    reinit!(v.cv, cell, x)
    for q_point in 1:getnquadpoints(v.cv)
        v.Δx[q_point] = spatial_coordinate(v.cv, q_point, x) - v.xc
    end
    return nothing
end
function Ferrite.shape_value(v::AffineCellValues{dim}, q_point::Int, i::Int) where {dim}
    nbase = getnbasefunctions(v.cv)
    i ≤ nbase && return shape_value(v.cv, q_point, i)
    k = i - nbase
    if k ≤ length(v.H_indices)
        (a, b) = v.H_indices[k]
        Δx_b = v.Δx[q_point][b]
        return Vec{dim}(l -> l == a ? Δx_b : 0.0)
    end
    c = v.g_indices[k - length(v.H_indices)]
    return Vec{dim}(l -> l == c ? 1.0 : 0.0)
end
function Ferrite.shape_gradient(v::AffineCellValues{dim}, q_point::Int, i::Int) where {dim}
    nbase = getnbasefunctions(v.cv)
    i ≤ nbase && return shape_gradient(v.cv, q_point, i)
    k = i - nbase
    k > length(v.H_indices) && return zero(Tensor{2, dim})
    (a, b) = v.H_indices[k]
    return Tensor{2, dim}((l, m) -> (l == a && m == b) ? 1.0 : 0.0)
end
function Ferrite.shape_hessian(v::AffineCellValues{dim}, q_point::Int, i::Int) where {dim}
    i ≤ getnbasefunctions(v.cv) && return shape_hessian(v.cv, q_point, i)
    return zero(Tensor{3, dim})
end

# Facet values with `n` additional base functions with zero value and gradient
struct PaddedFacetValues{FV <: FacetValues} <: Ferrite.AbstractFacetValues
    fv::FV
    n::Int
end
Ferrite.getnbasefunctions(v::PaddedFacetValues) = getnbasefunctions(v.fv) + v.n
Ferrite.getnquadpoints(v::PaddedFacetValues) = getnquadpoints(v.fv)
Ferrite.getdetJdV(v::PaddedFacetValues, q_point::Int) = getdetJdV(v.fv, q_point)
Ferrite.reinit!(v::PaddedFacetValues, cell, x::AbstractVector{<:Vec}, facetnr::Int) = reinit!(v.fv, cell, x, facetnr)
function Ferrite.shape_value(v::PaddedFacetValues, q_point::Int, i::Int)
    N = shape_value(v.fv, q_point, 1)
    return i ≤ getnbasefunctions(v.fv) ? shape_value(v.fv, q_point, i) : zero(N)
end

# Global map, T, from the dofs of `dh_attached`, where `H` (and `g`) are attached to `:u`, to the dofs of
# `dh_plain`, where u = ũ + H ⋅ (X - xc) + g in the nodes, such that a_plain = T * a_attached.
# Since linear Lagrange interpolations reproduce affine fields, the assembled attached system is
# Tᵀ * K_plain * T, and the residual Tᵀ * r_plain(T * a_attached).
function affine_map(dh_plain::DofHandler, dh_attached::DofHandler, xc::Vec{dim}, H_indices, g_indices = ()) where {dim}
    nplain = ndofs(dh_plain)
    @assert all(i -> i ≤ nplain, celldofs(dh_attached, 1)) # Spatial dofs are numbered first in both
    T = zeros(nplain, ndofs(dh_attached))
    for i in 1:nplain
        T[i, i] = 1.0
    end
    H_dofs = algebraic_dofs(dh_attached, :H)
    H_var = FerriteAssembly._algebraic_variable(dh_attached, :H)
    g_dofs = isempty(g_indices) ? Int[] : algebraic_dofs(dh_attached, :g)
    for sdh in dh_plain.subdofhandlers
        :u ∈ Ferrite.getfieldnames(sdh) || continue
        ru = dof_range(sdh, :u)
        for cellnr in sdh.cellset
            dofs = celldofs(dh_plain, cellnr)
            X = getcoordinates(dh_plain.grid, cellnr)
            for (k, i) in enumerate(ru)
                node, comp = divrem(k - 1, dim) .+ 1
                for (a, b) in H_indices
                    a == comp || continue
                    T[dofs[i], H_dofs[FerriteAssembly._algebraic_component(H_var, (a, b))]] = (X[node] - xc)[b]
                end
                for (m, a) in enumerate(g_indices)
                    a == comp && (T[dofs[i], g_dofs[a]] = 1.0)
                end
            end
        end
    end
    return T
end
