# Create FacetValues or CellValues automatically with minimal input 

"""
    autogenerate_facetvalues(fv::FacetValues, args...)

Just return the provided FacetValues 

    autogenerate_facetvalues(order::Int, ip_fun::Interpolation{RefShape}, ip_geo::Interpolation{RefShape})

Using quadrature rule, `fqr = FacetQuadratureRule{RefShape}(order)`, 
create `FacetValues(fqr, ip_fun, ip_geo)`
"""
autogenerate_facetvalues(fv::AbstractFacetValues, args...) = fv
function autogenerate_facetvalues(order::Int, ip::Interpolation{RefShape}, ip_geo::Interpolation{RefShape}) where RefShape
    return FacetValues(FacetQuadratureRule{RefShape}(order), ip, ip_geo)
end

"""
    autogenerate_cellvalues(cv::CellValues, args...)

Just return the provided CellValues 

    autogenerate_cellvalues(order::Int, ip_fun::Interpolation{RefShape}, ip_geo::Interpolation{RefShape})

Using quadrature rule, `qr = QuadratureRule{RefShape}(order)`, 
return `CellValues(qr, ip_fun, ip_geo)`
"""
autogenerate_cellvalues(cv::AbstractCellValues, args...) = cv
function autogenerate_cellvalues(order::Int, ip::Interpolation{RefShape}, ip_geo::Interpolation{RefShape}) where RefShape
    return CellValues(QuadratureRule{RefShape}(order), ip, ip_geo)
end
