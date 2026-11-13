"""
    FerriteAssembly.check_values_interface(values::AbstractCellValues, cell, x; hessians = false)

**Experimental**. Check that the cell values type of `values` supports the accessors that element
routines typically use, see [Values for attached variables](@ref). This is intended for testing
custom values types, e.g. for fields with attached algebraic variables (see
[`AttachAlgebraic`](@ref FerriteAssembly.AttachAlgebraic)). `values` are first reinitialized
with `reinit!(values, cell, x)`. Then, for all quadrature points and all base functions, it is
checked that
* `shape_value`, `shape_gradient`, and their derived accessors (`shape_symmetric_gradient`,
  `shape_divergence`, and `shape_curl` for vector-valued functions) can be called, and return the
  same type for all base functions;
* the corresponding `function_*` accessors can be called for local dof values `ae`, and
  for `ae` with a `dof_range`;
* `getdetJdV` and `spatial_coordinate(values, q_point, x)` can be called;
* none of these allocate.
If `hessians = true`, `shape_hessian` and `function_hessian` are included as well.

Returns `true` if all checks pass, and throws an error describing the first failure otherwise.
"""
function check_values_interface(values::Ferrite.AbstractCellValues, cell, x::AbstractVector{<:Vec}; hessians::Bool = false)
    reinit!(values, cell, x)
    n = getnbasefunctions(values)
    nqp = getnquadpoints(values)
    (n > 0 && nqp > 0) || error("The values must have at least one base function and one quadrature point")
    vector_valued = shape_value(values, 1, 1) isa Vec
    accessors = (shape_value, shape_gradient)
    if vector_valued
        accessors = (accessors..., shape_symmetric_gradient, shape_divergence, shape_curl)
    end
    if hessians
        accessors = (accessors..., shape_hessian)
    end
    for accessor in accessors
        T = typeof(accessor(values, 1, 1))
        for q_point in 1:nqp, i in 1:n
            val = _call_accessor(accessor, values, q_point, i)
            if typeof(val) !== T
                error("$accessor(values, $q_point, $i) returns a $(typeof(val)), but $accessor(values, 1, 1) returns a $T")
            end
        end
    end
    ae = collect(range(0.1, 1.0; length = n))
    ae_padded = vcat(0.0, ae, 0.0)
    flags = (Val(vector_valued), Val(hessians))
    for q_point in 1:nqp
        r1 = _call_functions(values, q_point, ae, eachindex(ae), flags...)
        r2 = _call_functions(values, q_point, ae_padded, 2:(n + 1), flags...)
        r1 == r2 || error("The function values differ when `ae` is given with and without a `dof_range`")
        getdetJdV(values, q_point)
        spatial_coordinate(values, q_point, x)
    end
    _accessor_loop(values, x, ae, flags...) # Compile
    nalloc = @allocated _accessor_loop(values, x, ae, flags...)
    nalloc == 0 || error("Calling the accessors allocates ($nalloc bytes)")
    return true
end

function _call_accessor(accessor::F, values, q_point, i) where {F}
    try
        return accessor(values, q_point, i)
    catch
        @error "$accessor(values, $q_point, $i) failed"
        rethrow()
    end
end

function _call_functions(values, q_point, ae, dofrange, ::Val{vector_valued}, ::Val{hessians}) where {vector_valued, hessians}
    vals = (function_value(values, q_point, ae, dofrange), function_gradient(values, q_point, ae, dofrange))
    vector_vals = vector_valued ? (
        function_symmetric_gradient(values, q_point, ae, dofrange),
        function_divergence(values, q_point, ae, dofrange),
        function_curl(values, q_point, ae, dofrange)) : ()
    hessian_vals = hessians ? (Ferrite.function_hessian(values, q_point, ae, dofrange),) : ()
    return (vals..., vector_vals..., hessian_vals...)
end

function _accessor_loop(values, x, ae, vv::Val{vector_valued}, hv::Val{hessians}) where {vector_valued, hessians}
    s = 0.0
    for q_point in 1:getnquadpoints(values)
        s += getdetJdV(values, q_point)
        s += norm(spatial_coordinate(values, q_point, x))
        for i in 1:getnbasefunctions(values)
            s += norm(shape_value(values, q_point, i)) + norm(shape_gradient(values, q_point, i))
            if vector_valued
                s += norm(shape_symmetric_gradient(values, q_point, i)) + shape_divergence(values, q_point, i)
                s += norm(shape_curl(values, q_point, i))
            end
            hessians && (s += norm(shape_hessian(values, q_point, i)))
        end
        s += sum(map(norm, _call_functions(values, q_point, ae, eachindex(ae), vv, hv)))
    end
    return s
end
