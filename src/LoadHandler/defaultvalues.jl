# Create FaceValues or CellValues automatically with minimial input 

"""
    autogenerate_facevalues(fv::FaceValues, args...)

Just return the provided FaceValues 

    autogenerate_facevalues(order::Int, ip_fun::Interpolation{RefShape}, ip_geo::Interpolation{RefShape})

Using quadrature rule, `fqr = FaceQuadratureRule{RefShape}(order)`, 
create `FaceValues(fqr, ip_fun, ip_geo)`
"""
autogenerate_facevalues(fv::FaceValues, args...) = fv
function autogenerate_facevalues(order::Int, ip::Interpolation{RefShape}, ip_geo::Interpolation{RefShape}) where RefShape
    return FaceValues(FaceQuadratureRule{RefShape}(order), ip, ip_geo)
end

"""
    autogenerate_cellvalues(cv::CellValues, args...)

Just return the provided CellValues 

    autogenerate_cellvalues(order::Int, ip_fun::Interpolation{RefShape}, ip_geo::Interpolation{RefShape})

Using quadrature rule, `qr = QuadratureRule{RefShape}(order)`, 
return `CellValues(qr, ip_fun, ip_geo)`
"""
autogenerate_cellvalues(cv::CellValues, args...) = cv
function autogenerate_cellvalues(order::Int, ip::Interpolation{RefShape}, ip_geo::Interpolation{RefShape}) where RefShape
    return CellValues(QuadratureRule{RefShape}(order), ip, ip_geo)
end
