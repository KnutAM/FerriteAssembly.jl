"""
    BodyLoad(field_name::Symbol, cv_info::Union{CellValues,QuadratureRule,Int}, cellset::Set{Int}, f)
    BodyLoad(field_name::Symbol, cv_info::Union{CellValues,QuadratureRule,Int}, f)

Define a body load contribution with the weak forms according to 
```math
\\int_{\\Omega} f \\ \\delta u \\ \\mathrm{d}\\Omega, \\quad \\text{Scalar fields} \\\\

\\int_{\\Omega} \\boldsymbol{f} \\cdot \\boldsymbol{\\delta u} \\ \\mathrm{d}\\Omega,
\\quad \\text{Vector fields} 
```
where ``\\Omega`` is the region where the contribution is active. 
``f``, or ``\\boldsymbol{f}``, is the prescribed volumetric contribution (e.g. force/volume), 
defined by a function with signatures

`f(x::Vec, t::Real) -> Number` (Scalar field)

`f(x::Vec, t::Real) -> Vec{dim}` (Vector field)

where `x` is the spatial position of the current quadrature point and `t` is the 
current time. The remaining input arguments are

* `fieldname` describes the field on which the load should act
* `cv_info` gives required input to determine the cellvalues. The following input types are accepted:
  - `Int` giving the integration order to use. `CellValues` are deduced from the interpolation 
    of `fieldname` and the output of `f`.   
  - `QuadratureRule` matching the interpolation for `fieldname` for cells in `cellset`. 
    `CellValues` are deduced from the output of `f`
  - `CellValues` matching the interpolation for `fieldname` for the cells in `cellset` and output of `f`

* `cellset` describes which cells the load is applied to. If not given, the load is applied to all cells 
"""
struct BodyLoad{CVI,FUN}
    fieldname::Symbol
    cv_info::CVI
    cellset::Union{Set{Int},Nothing}
    f::FUN # f(x::Vec, time)->{cv::CellScalarValues ? Number : Vec}
end
BodyLoad(fieldname::Symbol, cv_info, f) = BodyLoad(fieldname, cv_info, nothing, f)

struct BodyLoadMaterial{FUN}
    f::FUN
    dr::UnitRange{Int}
end

function element_residual!(fe::Vector, ::Any, ::Vector, m::BodyLoadMaterial, cv::CellValues, cellbuffer)
    checkbounds(fe, m.dr)
    t = get_time_increment(cellbuffer) # Abuse...
    for q_point in 1:getnquadpoints(cv)
        dΩ = getdetJdV(cv, q_point)
        x = spatial_coordinate(cv, q_point, getcoordinates(cellbuffer))
        b = m.f(x, t)
        for (i, I) in pairs(m.dr)
            δu = shape_value(cv, q_point, i)
            @inbounds fe[I] += (δu ⋅ b) * dΩ
        end
    end
end

function add_bodyload!(bodyloads::Dict{String}, bodyload::BodyLoad, dh::DofHandler)
    add_bodyload!(bodyloads, bodyload, SubDofHandler(dh))
end

function add_bodyload!(bodyloads::Dict{String,BT}, bodyload::BodyLoad, sdh::SubDofHandler) where BT
    material = BodyLoadMaterial(bodyload.f, dof_range(sdh, bodyload.fieldname))

    ip = Ferrite.getfieldinterpolation(sdh, bodyload.fieldname)
    ip_geo = Ferrite.default_interpolation(getcelltype(sdh))
    cv = autogenerate_cellvalues(bodyload.cv_info, ip, ip_geo, bodyload.f)

    set = bodyload.cellset===nothing ? getcellset(sdh) : bodyload.cellset
    
    domain_spec = DomainSpec(sdh, material, cv; set=set)
    threading = BT <: ThreadedDomainBuffer
    bodyloads[string(length(bodyloads)+1)] = setup_domainbuffer(domain_spec; threading=threading)
end

function add_bodyload!(bodyloads::Dict{String}, bodyload::BodyLoad, dh::MixedDofHandler)
    contribution = false
    for fh in dh.fieldhandlers
        overlaps = overlaps_with_cellset(bodyload.cellset, fh.cellset)
        if overlaps && bodyload.fieldname ∈ Ferrite.getfieldnames(fh)
            contribution = true
            add_bodyload!(bodyloads, bodyload, SubDofHandler(dh, fh))
        end
    end
    contribution || @warn "No contributions added to the LoadHandler"
end
