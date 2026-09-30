# Create FacetValues or CellValues automatically with minimal input 

"""
    autogenerate_facetvalues(fv::AbstractFacetValues, args...)

Just return the provided FacetValues

    autogenerate_facetvalues(order::Int, ip_fun::Interpolation{RefShape}, ip_geo::Interpolation{RefShape})

Using quadrature rule, `fqr = FacetQuadratureRule{RefShape}(order)`,
create `FacetValues(fqr, ip_fun, ip_geo)`

    autogenerate_facetvalues(fqr::FacetQuadratureRule{RefShape}, ip_fun::Interpolation{RefShape}, ip_geo::Interpolation{RefShape})

Create `FacetValues(fqr, ip_fun, ip_geo)` directly from the given quadrature rule.
"""
autogenerate_facetvalues(fv::AbstractFacetValues, args...) = fv
function autogenerate_facetvalues(order::Int, ip::Interpolation{RefShape}, ip_geo::Interpolation{RefShape}) where RefShape
    return FacetValues(FacetQuadratureRule{RefShape}(order), ip, ip_geo)
end
function autogenerate_facetvalues(fqr::FacetQuadratureRule{RefShape}, ip::Interpolation{RefShape}, ip_geo::Interpolation{RefShape}) where RefShape
    return FacetValues(fqr, ip, ip_geo)
end

"""
    autogenerate_cellvalues(cv::AbstractCellValues, args...)

Just return the provided CellValues

    autogenerate_cellvalues(order::Int, ip_fun::Interpolation{RefShape}, ip_geo::Interpolation{RefShape})

Using quadrature rule, `qr = QuadratureRule{RefShape}(order)`,
return `CellValues(qr, ip_fun, ip_geo)`

    autogenerate_cellvalues(qr::QuadratureRule{RefShape}, ip_fun::Interpolation{RefShape}, ip_geo::Interpolation{RefShape})

Create `CellValues(qr, ip_fun, ip_geo)` directly from the given quadrature rule.
"""
autogenerate_cellvalues(cv::AbstractCellValues, args...) = cv
function autogenerate_cellvalues(order::Int, ip::Interpolation{RefShape}, ip_geo::Interpolation{RefShape}) where RefShape
    return CellValues(QuadratureRule{RefShape}(order), ip, ip_geo)
end
function autogenerate_cellvalues(qr::QuadratureRule{RefShape}, ip::Interpolation{RefShape}, ip_geo::Interpolation{RefShape}) where RefShape
    return CellValues(qr, ip, ip_geo)
end
