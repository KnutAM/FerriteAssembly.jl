"""
    Neumann(field_name::Symbol, fv_info::Union{FaceValues,QuadratureRule,Int}, faceset::Set{FaceIndex}, f)

Define a Neumann contribution with the weak forms according to 
```math
\\int_{\\Gamma} f \\ \\delta u \\ \\mathrm{d}\\Gamma, \\quad \\text{Scalar fields} \\\\

\\int_{\\Gamma} \\boldsymbol{f} \\cdot \\boldsymbol{\\delta u} \\ \\mathrm{d}\\Gamma,
\\quad \\text{Vector fields} 
```
where ``\\Gamma`` is the boundary where the contribution is active. 
``f``, or ``\\boldsymbol{f}``, is the prescribed Neumann value, 
defined by a function with signatures

`f(x::Vec, time, n::Vec) -> Number` (Scalar field)

`f(x::Vec, time, n::Vec) -> Vec{dim}` (Vector field)

where `x` is the spatial position of the current quadrature point, `time` is the 
current time, and `n` is the face normal vector. The remaining input arguments are

* `fieldname` describes the field on which the boundary condition should abstract
* `fv_info` gives required input to determine the facevalues. The following input types are accepted:
  - `FaceValues` matching the interpolation for `fieldname` for the faces in `faceset` and output of `f`
  - `QuadratureRule` matching the interpolation for `fieldname` for faces in `faceset` `FaceValues` are deduced 
    from the output of `f`
  - `Int` giving the integration order to use. `FaceValues` are deduced from the interpolation 
    of `fieldname` and the output of `f`. 
* `faceset` describes which faces the BC is applied to
"""
struct Neumann{FVI,FUN}
    fieldname::Symbol
    fv_info::FVI
    faceset::Set{FaceIndex}
    f::FUN # f(x::Vec, time, n::Vec)->{FV::FaceScalarValues ? Number : Vec}
end

# Internal
struct NeumannMaterial{FUN}
    f::FUN
    dr::UnitRange{Int}
end
function face_residual!(fe::Vector, ::Vector, m::NeumannMaterial, fv::FaceValues, facebuffer)
    checkbounds(fe, m.dr)
    t = get_time_increment(facebuffer) # Abuse...
    for q_point in 1:getnquadpoints(fv)
        dΓ = getdetJdV(fv, q_point)
        x = spatial_coordinate(fv, q_point, getcoordinates(facebuffer))
        n = getnormal(fv, q_point)
        b = m.f(x, t, n)
        for (i, I) in pairs(m.dr)
            δu = shape_value(fv, q_point, i)
            @inbounds fe[I] += (δu ⋅ b) * dΓ
        end
    end
end

function add_neumann!(nbcs::Dict{String}, nbc::Neumann, dh::DofHandler)
    return add_neumann!(nbcs, nbc, SubDofHandler(dh))
end

function add_neumann!(nbcs::Dict{String,BT}, nbc::Neumann, sdh::SubDofHandler) where BT
    material = NeumannMaterial(nbc.f, dof_range(sdh, nbc.fieldname))

    ip = Ferrite.getfieldinterpolation(sdh, nbc.fieldname)
    ip_geo = Ferrite.default_interpolation(getcelltype(sdh))
    fv = get_facevalues(nbc.fv_info, ip, ip_geo, nbc.f)
    
    domain_spec = DomainSpec(sdh, material, fv; set=nbc.faceset)
    threading = BT <: ThreadedDomainBuffer
    nbcs[string(length(nbcs)+1)] = setup_domainbuffer(domain_spec; threading=threading)
end

function add_neumann!(nbcs::Dict{String}, nbc::Neumann, dh::MixedDofHandler)
    contribution = false
    for fh in dh.fieldhandlers
        overlaps = overlaps_with_cellset(nbc.faceset, fh.cellset)
        if overlaps && nbc.fieldname ∈ Ferrite.getfieldnames(fh)
            contribution = true
            add_neumann!(nbcs, nbc, SubDofHandler(dh, fh))
        end
    end
    contribution || @warn "No contributions added to the LoadHandler"
end
