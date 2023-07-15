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

`f(x::Vec, time) -> Number` (Scalar field)

`f(x::Vec, time) -> Vec{dim}` (Vector field)

where `x` is the spatial position of the current quadrature point and `time` is the 
current time. The remaining input arguments are

* `fieldname` describes the field on which the load should act
* `cv_info` gives required input to determine the cellvalues. The following input types are accepted:
  - `CellValues` matching the interpolation for `fieldname` for the cells in `cellset` and output of `f`
  - `QuadratureRule` matching the interpolation for `fieldname` for cells in `cellset`. 
    `CellValues` are deduced from the output of `f`
  - `Int` giving the integration order to use. `CellValues` are deduced from the interpolation 
    of `fieldname` and the output of `f`. 
* `cellset` describes which cells the load is applied to. If not given, the load is applied to all cells 
"""
struct BodyLoad{CVI,FUN}
    fieldname::Symbol
    cv_info::CVI
    cellset::Union{Set{Int},Nothing}
    f::FUN # f(x::Vec, time)->{cv::CellScalarValues ? Number : Vec}
end
BodyLoad(fieldname::Symbol, cv_info, f) = BodyLoad(fieldname, cv_info, nothing, f)

# Internal
struct BodyLoadData{CV,FUN}
    fieldname::Symbol   # Only for information 
    dofrange::UnitRange{Int}
    cellvalues::CV
    cellset::Set{Int}
    f::FUN
end

function BodyLoadData(dh::DofHandler, spec::BodyLoad)
    cellset = isnothing(spec.cellset) ? Set(1:getncells(dh.grid)) : spec.cellset
    cell = getcells(dh.grid, first(cellset))
    BodyLoadData(dh, spec, cellset, cell)
end

function BodyLoadData(dh_fh::Union{DofHandler,FieldHandler}, spec::BodyLoad, cellset::Set{Int}, ::C) where C<:Ferrite.AbstractCell
    dofrange = dof_range(dh_fh, spec.fieldname)
    ip = Ferrite.getfieldinterpolation(dh_fh, Ferrite.find_field(dh_fh, spec.fieldname))
    ip_geo = Ferrite.default_interpolation(C)
    cv = get_cellvalues(spec.cv_info, ip, ip_geo, spec.f)
    return BodyLoadData(spec.fieldname, dofrange, cv, cellset, spec.f)
end

function add_bodyload!(bodyloads::Vector, bodyload::BodyLoad, dh::DofHandler)
    push!(bodyloads, BodyLoadData(dh, bodyload))
end

_intersect_nothing(::Nothing, set) = set 
_intersect_nothing(a, b) = intersect(a, b)

function add_bodyload!(bodyloads::Vector, bodyload::BodyLoad, dh::MixedDofHandler)
    contribution = false
    for fh in dh.fieldhandlers
        cellset = _intersect_nothing(bodyload.cellset, fh.cellset)
        if !isempty(cellset) && bodyload.fieldname ∈ Ferrite.getfieldnames(fh)
            contribution = true
            cell = getcells(dh.grid, first(cellset))
            push!(bodyloads, BodyLoadData(fh, bodyload, cellset, cell))
        end
    end
    contribution || @warn "No contributions added to the NeumannHandler"
end

function apply_bodyload!(f::Vector{T}, bodyload::BodyLoadData, dh, time) where T
    dofs = collect(bodyload.dofrange)
    fe = zeros(T, length(dofs))
    for cell in CellIterator(dh, bodyload.cellset)
        calculate_bodyload_contribution!(fe, cell, bodyload.cellvalues, time, bodyload.f)
        assemble!(f, view(celldofs(cell), dofs), fe)
    end
end

function calculate_bodyload_contribution!(fe::Vector, cell::CellCache, cv::CellValues, time, f)
    fill!(fe, 0)
    reinit!(cv, cell)
    for q_point in 1:getnquadpoints(cv)
        dΩ = getdetJdV(cv, q_point)
        x = spatial_coordinate(cv, q_point, getcoordinates(cell))
        b = f(x, time)
        for i in 1:getnbasefunctions(cv)
            δu = shape_value(cv, q_point, i)
            fe[i] += (δu ⋅ b) * dΩ
        end
    end
end
