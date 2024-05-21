"""
    Neumann(field_name::Symbol, fv_info::Union{FacetValues,QuadratureRule,Int}, facetset::AbstractSet{FacetIndex}, f)

Define a Neumann contribution with the weak forms according to 
```math
\\int_{\\Gamma} f \\ \\delta u \\ \\mathrm{d}\\Gamma, \\quad \\text{Scalar fields} \\\\

\\int_{\\Gamma} \\boldsymbol{f} \\cdot \\boldsymbol{\\delta u} \\ \\mathrm{d}\\Gamma,
\\quad \\text{Vector fields} 
```
where ``\\Gamma`` is the boundary where the contribution is active. 
``f``, or ``\\boldsymbol{f}``, is the prescribed Neumann value, 
defined by a function with signatures

`f(x::Vec, t::Real, n::Vec) -> Number` (Scalar field)

`f(x::Vec, t::Real, n::Vec) -> Vec{dim}` (Vector field)

where `x` is the spatial position of the current quadrature point, `t` is the 
current time, and `n` is the facet normal vector. The remaining input arguments are

* `fieldname` describes the field on which the boundary condition should abstract
* `fv_info` gives required input to determine the facetvalues. The following input types are accepted:
  - `Int` giving the integration order to use. `FacetValues` are deduced from the interpolation 
    of `fieldname` and the output of `f`. 
  - `QuadratureRule` matching the interpolation for `fieldname` for facets in `facetset` `FacetValues` are deduced 
    from the output of `f`
  - `FacetValues` matching the interpolation for `fieldname` for the facets in `facetset` and output of `f`
* `facetset` describes which facets the BC is applied to
"""
struct Neumann{FVI,FUN}
    fieldname::Symbol
    fv_info::FVI
    facetset::AbstractSet{FacetIndex}
    f::FUN # f(x::Vec, time, n::Vec)->{ip::ScalarInterpolation ? Number : Vec}
end

# Internal
struct NeumannMaterial{FUN}
    f::FUN
    dr::UnitRange{Int}
end
function facet_residual!(fe::Vector, ::Vector, m::NeumannMaterial, fv::FacetValues, facetbuffer)
    checkbounds(fe, m.dr)
    t = get_time_increment(facetbuffer) # Abuse...
    for q_point in 1:getnquadpoints(fv)
        dΓ = getdetJdV(fv, q_point)
        x = spatial_coordinate(fv, q_point, getcoordinates(facetbuffer))
        n = getnormal(fv, q_point)
        b = m.f(x, t, n)
        for (i, I) in pairs(m.dr)
            δu = shape_value(fv, q_point, i)
            @inbounds fe[I] += (δu ⋅ b) * dΓ
        end
    end
end

function add_neumann!(nbcs::Dict{String,BT}, nbc::Neumann, sdh::SubDofHandler) where BT
    material = NeumannMaterial(nbc.f, dof_range(sdh, nbc.fieldname))

    ip = Ferrite.getfieldinterpolation(sdh, nbc.fieldname)
    ip_geo = Ferrite.default_interpolation(getcelltype(sdh))
    fv = autogenerate_facetvalues(nbc.fv_info, ip, ip_geo)
    
    domain_spec = DomainSpec(sdh, material, fv; set=nbc.facetset)
    threading = BT <: ThreadedDomainBuffer
    nbcs[string(length(nbcs)+1)] = setup_domainbuffer(domain_spec; threading=threading)
end

function add_neumann!(nbcs::Dict{String}, nbc::Neumann, dh::DofHandler)
    contribution = false
    for sdh in dh.subdofhandlers
        overlaps = overlaps_with_cellset(nbc.facetset, _getcellset(sdh))
        if overlaps && nbc.fieldname ∈ Ferrite.getfieldnames(sdh)
            contribution = true
            add_neumann!(nbcs, nbc, sdh)
        end
    end
    contribution || @warn "No contributions added to the LoadHandler"
end
