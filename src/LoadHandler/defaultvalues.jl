# Create FaceValues or CellValues automatically with minimial input 

"""
    autogenerate_facevalues(fv::FaceValues, args...)

Just return the provided FaceValues 

    autogenerate_facevalues(qr::QuadratureRule, ip_fun::Interpolation, ip_geo::Interpolation, f)

Depending on the return type of `f(x,t,n)`, return `FaceScalarValues(qr, ip_fun, ip_geo)`
or `FaceVectorValues(qr, ip_fun, ip_geo)`.

    autogenerate_facevalues(order::Int, ip_fun::Interpolation{dim,RefShape}, ip_geo::Interpolation{dim,RefShape}, f)

Using quadrature rule, `qr = QuadratureRule{dim-1,RefShape}(order)`, call `autogenerate_facevalues(qr, ip_fun, ip_geo, f)`
"""
function autogenerate_facevalues end

# If a facevalue has already been given, use this value 
autogenerate_facevalues(fv::FaceValues, args...) = fv

# Create default quadrule of given order
function autogenerate_facevalues(order::Int, ip::Interpolation{dim,RefShape}, ip_geo::Interpolation{dim,RefShape}, f) where {dim, RefShape}
    return autogenerate_facevalues(QuadratureRule{dim-1,RefShape}(order), ip, ip_geo, f)
end

# Use the given function to determine if the output should be a scalar or vector. 
function autogenerate_facevalues(qr::QuadratureRule, ip::Interpolation{dim}, ip_geo, f) where dim
    return autogenerate_facevalues(qr, ip, ip_geo, f(zero(Vec{dim}), 0.0, zero(Vec{dim})))
end
function autogenerate_facevalues(qr::QuadratureRule{<:Any,RefShape}, ip::Interpolation{<:Any,RefShape}, 
                        ip_geo::Interpolation{<:Any,RefShape}, ::Vec) where RefShape
    return FaceVectorValues(qr, ip, ip_geo)
end
function autogenerate_facevalues(qr::QuadratureRule{<:Any,RefShape}, ip::Interpolation{<:Any,RefShape}, 
                        ip_geo::Interpolation{<:Any,RefShape}, ::Number) where RefShape
    return FaceScalarValues(qr, ip, ip_geo)
end
#= Not required, will throw assertion error upon FaceValues construction (but keep if later decided a nicer error is good)
function autogenerate_facevalues(qr::QuadratureRule, ip::Interpolation, ip_geo::Interpolation, fval::Union{Vec,Number})
    throw(ArgumentError("qr, $(typeof(qr)), ip, $(typeof(ip)), and ip_geo, $(typeof(ip_geo)), doesn't seem compatible. (info: fval=$fval)"))
end
=#

"""
    autogenerate_cellvalues(fv::CellValues, args...)

Just return the provided CellValues 

    autogenerate_cellvalues(qr::QuadratureRule, ip_fun::Interpolation, ip_geo::Interpolation, f)

Depending on the return type of `f(x,t)`, return `CellScalarValues(qr, ip_fun, ip_geo)`
or `CellVectorValues(qr, ip_fun, ip_geo)`.

    autogenerate_cellvalues(order::Int, ip_fun::Interpolation{dim,RefShape}, ip_geo::Interpolation{dim,RefShape}, f)

Using quadrature rule, `qr = QuadratureRule{dim,RefShape}(order)`, call `autogenerate_cellvalues(qr, ip_fun, ip_geo, f)`
"""
function autogenerate_cellvalues end

# If a cellvalue has already been given, use this value 
autogenerate_cellvalues(cv::CellValues, args...) = cv

# Create default quadrule of given order
function autogenerate_cellvalues(order::Int, ip::Interpolation{dim,RefShape}, ip_geo::Interpolation{dim,RefShape}, f) where {dim, RefShape}
    return autogenerate_cellvalues(QuadratureRule{dim,RefShape}(order), ip, ip_geo, f)
end

# Use the given function to determine if the output should be a scalar or vector. 
function autogenerate_cellvalues(qr::QuadratureRule, ip::Interpolation{dim}, ip_geo, f) where dim
    return autogenerate_cellvalues(qr, ip, ip_geo, f(zero(Vec{dim}), 0.0))
end
function autogenerate_cellvalues(qr::QuadratureRule{dim,RefShape}, ip::Interpolation{dim,RefShape}, 
                        ip_geo::Interpolation{dim,RefShape}, ::Vec) where {dim, RefShape}
    return CellVectorValues(qr, ip, ip_geo)
end
function autogenerate_cellvalues(qr::QuadratureRule{dim,RefShape}, ip::Interpolation{dim,RefShape}, 
                        ip_geo::Interpolation{dim,RefShape}, ::Number) where {dim, RefShape}
    return CellScalarValues(qr, ip, ip_geo)
end
#= Not required, will throw assertion error upon CellValues construction (but keep if later decided a nicer error is good)
function autogenerate_cellvalues(qr::QuadratureRule, ip::Interpolation, ip_geo::Interpolation, fval::Union{Vec,Number})
    throw(ArgumentError("qr, $(typeof(qr)), ip, $(typeof(ip)), and ip_geo, $(typeof(ip_geo)), doesn't seem compatible. (info: fval=$fval)"))
end
=#