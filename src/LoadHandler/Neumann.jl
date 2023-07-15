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
struct NeumannData{FV,FUN}
    fieldname::Symbol   # Only for information 
    dofrange::UnitRange{Int}
    facevalues::FV 
    faceset::Set{FaceIndex}
    f::FUN
end

function NeumannData(dh::DofHandler, spec::Neumann)
    cell = getcells(dh.grid, first(first(spec.faceset)))
    NeumannData(dh, spec, spec.faceset, cell)
end

function NeumannData(dh_fh::Union{DofHandler,FieldHandler}, spec::Neumann, faceset::Set{FaceIndex}, ::C) where C<:Ferrite.AbstractCell
    dofrange = dof_range(dh_fh, spec.fieldname)
    ip = Ferrite.getfieldinterpolation(dh_fh, Ferrite.find_field(dh_fh, spec.fieldname))
    ip_geo = Ferrite.default_interpolation(C)
    fv = get_facevalues(spec.fv_info, ip, ip_geo, spec.f)
    return NeumannData(spec.fieldname, dofrange, fv, faceset, spec.f)
end

function add_neumann!(nbcs::Vector, nbc::Neumann, dh::DofHandler)
    push!(nbcs, NeumannData(dh, nbc))
end

function add_neumann!(nbcs::Vector, nbc::Neumann, dh::MixedDofHandler)
    contribution = false
    for fh in dh.fieldhandlers
        faceset = intersect_with_cellset(nbc.faceset, fh.cellset)
        if !isempty(faceset) && nbc.fieldname ∈ Ferrite.getfieldnames(fh)
            contribution = true
            cell = getcells(dh.grid, first(first(faceset)))
            push!(nbcs, NeumannData(fh, nbc, faceset, cell))
        end
    end
    contribution || @warn "No contributions added to the NeumannHandler"
end

function intersect_with_cellset(faceset::Set{FaceIndex}, cellset)
    return Set(face for face in faceset if first(face) in cellset)
end

function apply_neumann!(f::Vector{T}, nbc::NeumannData, dh, time) where T
    dofs = collect(nbc.dofrange)
    fe = zeros(T, length(dofs))
    for face in FaceIterator(dh, nbc.faceset)
        calculate_neumann_contribution!(fe, face, nbc.facevalues, time, nbc.f)
        assemble!(f, view(celldofs(face), dofs), fe)
    end
end

function calculate_neumann_contribution!(fe::Vector, face::FaceCache, fv::FaceValues, time, f)
    fill!(fe, 0)
    reinit!(fv, face)
    for q_point in 1:getnquadpoints(fv)
        dΓ = getdetJdV(fv, q_point)
        x = spatial_coordinate(fv, q_point, getcoordinates(face))
        n = getnormal(fv, q_point)
        b = f(x, time, n)
        for i in 1:getnbasefunctions(fv)
            δu = shape_value(fv, q_point, i)
            fe[i] += (δu ⋅ b) * dΓ
        end
    end
end
