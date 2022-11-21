abstract type AbstractCellBuffer end

struct CellBuffer{dim,T,CV,MT,CL,CT} <: AbstractCellBuffer
    ae_old::Vector{T}           # Old element dof values
    ae::Vector{T}               # New element dof values
    re::Vector{T}               # Residual/force vector 
    Ke::Matrix{T}               # Element stiffness matrix
    dofs::Vector{Int}           # celldofs
    coords::Vector{Vec{dim,T}}  # cellcoords
    cellvalues::CV              # Can also be e.g. tuple of cellvalues
    material::MT                # User material definition
    cell_load::CL               # User source term/body force definition
    cache::CT                   # CellBuffer as user pleases
end

"""
    CellBuffer(
        numdofs::Int, numnodes::Int, ::Val{dim}, 
        cellvalues, material, cell_load=nothing, cache=nothing) -> CellBuffer

Create a cell cache for an element with `numdofs` degrees of freedom and
`numnodes` nodes with dimension `dim`. Add the given `cellvalues`, `material`, 
and `cache` to the `CellBuffer` as well. 

    CellBuffer(
        dh::DofHandler, 
        cellvalues, material, cell_load=nothing, cache=nothing) -> CellBuffer

Use `dh` to get `numdofs`, `numnodes`, and `dim`, before calling the above method definition. 

    CellBuffer(
        dh::MixedDofHandler, fh::FieldHandler, 
        cellvalues, material, cell_load=nothing, cache=nothing) -> CellBuffer

Use `dh` and `fh` to get `numdofs`, `numnodes`, and `dim`, 
before calling the first `CellBuffer` method definition. 

    CellBuffer(
        dh::MixedDofHandler, cellvalues, 
        material, cell_load=nothing, cache=nothing) -> NTuple{N,<:CellBuffer}
    
Return a tuple of `CellBuffer`s for each `FieldHandler` in `dh.fieldhandlers`.
`cellvalues[i]` corresponds to `dh.fieldhandlers[i]`, and so does 
`materials[i]`, `cell_load[i]` and `caches[i]`. If only one `CellValues`, `material`, `cell_load`, and/or `cache`
is given (not as a `::Tuple`), the same is used for all `fieldhandlers`. 
If a tuple of `cellvalues` (or materials/cell_load/caches) should be used for each cell, 
and the same tuple should be used for each fieldhandler, 
then it must be given as a tuple of tuples. 
(Often, it is better to give a `NamedTuple` of e.g. `CellValues` to be used for every fieldhandler)

The `Ferrite.jl` functions `getcoordinates(::CellBuffer)` and `celldofs(::CellBuffer)` are defined and can 
be used inside an element routine to get the current cell's coordinates and dof-numbers. 
"""
function CellBuffer(numdofs::Int, numnodes::Int, ::Val{dim}, cellvalues, material, cell_load=nothing, cache=nothing) where dim
    return CellBuffer(
        zeros(numdofs), zeros(numdofs), zeros(numdofs), zeros(numdofs,numdofs), 
        zeros(Int, numdofs), zeros(Vec{dim}, numnodes), 
        cellvalues, material, cell_load, cache
        )
end

function CellBuffer(dh::DofHandler{dim}, args...) where dim
    return CellBuffer(ndofs_per_cell(dh), Ferrite.nnodes_per_cell(dh), Val{dim}(), args...)
end

function CellBuffer(dh::MixedDofHandler{dim}, fh::FieldHandler, args...) where dim
    return CellBuffer(ndofs_per_cell(dh, fh), Ferrite.nnodes_per_cell(dh, fh), Val{dim}(), args...)
end

function CellBuffer(dh::MixedDofHandler, cellvalues, materials, cell_load=nothing, caches=nothing)
    numfh = length(dh.fieldhandlers)
    cellvalues_ = _maketuple(cellvalues, numfh)
    materials_ = _maketuple(materials, numfh)
    cell_load_ = _maketuple(cell_load, numfh)
    caches_ = _maketuple(caches, numfh)
    return ntuple(i->CellBuffer(dh, dh.fieldhandlers[i], cellvalues_[i], materials_[i], cell_load_[i], caches_[i]), numfh)
end

# Need special treatment due to type ambiguity
CellBuffer(dh::DofHandler, cv, mtrls::Dict{String}, args...; kwargs...) = _dictCellBuffer(dh, cv, mtrls, args...; kwargs...)
CellBuffer(dh::MixedDofHandler, cv, mtrls::Dict{String}, args...; kwargs...) = _dictCellBuffer(dh, cv, mtrls, args...; kwargs...)
function _dictCellBuffer(dh::Union{DofHandler,MixedDofHandler}, cv, materials::Dict{String}, cell_load=nothing, caches=nothing)
    setnames = keys(materials)
    cv_ = _makedict(cv, setnames)
    cell_load_ = _makedict(cell_load, setnames)
    caches_ = _makedict(caches, setnames)
    return Dict(setname=>
        CellBuffer(dh, cv_[setname], materials[setname], cell_load_[setname], caches_[setname])
        for setname in setnames)
end


# Convenience access functions to CellBuffer for use inside element routines
# Currently only Ferrite overloads used, because I currently feel that the rest are only 
# convenient if exported, and in that case it clotters namespace more than the convenience motivates

@inline Ferrite.getcoordinates(c::CellBuffer) = c.coords
#=
@inline get_ae(c::CellBuffer) = c.ae
@inline get_aeold(c::CellBuffer) = c.ae_old
@inline get_load(c::CellBuffer) = c.cell_load
@inline get_cache(c::CellBuffer) = c.cache

get_ae_view(c, dh_fh::Union{DofHandler,FieldHandler}, field::Symbol) = view(get_ae(c), dof_range(dh_fh, field))
get_aeold_view(c, dh_fh::Union{DofHandler,FieldHandler}, field::Symbol) = view(get_aeold(c), dof_range(dh_fh, field))
get_re_view(c, dh_fh::Union{DofHandler,FieldHandler}, field::Symbol) = view(get_re(c), dof_range(dh_fh, field))
=#

# Required functions for a custom CellBuffer
@inline get_Ke(c::CellBuffer) = c.Ke
@inline get_re(c::CellBuffer) = c.re
@inline get_ae(c::CellBuffer) = c.ae
@inline get_material(c::CellBuffer) = c.material
@inline get_cellvalues(c::CellBuffer) = c.cellvalues
@inline Ferrite.celldofs(c::CellBuffer) = c.dofs

"""
    Ferrite.reinit!(c::CellBuffer, dh::AbstractDofHandler, cellnum::Int, anew, aold)

Reinitialize the `c::CellBuffer` for cell number `cellnum`.
The global degree of freedom vectors `anew` (current) and `aold` are used
to update the cell degree of freedom vectors in `c`.
If the global vectors are instead `::Nothing`, the corresponding cell values are set to `NaN`
The element stiffness, `c.Ke`, and residual, `c.re`, are also zeroed. 
"""
function Ferrite.reinit!(c::CellBuffer, dh::Ferrite.AbstractDofHandler, cellnum::Int, anew, aold)
    celldofs!(c.dofs, dh, cellnum)
    Ferrite.cellcoords!(c.coords, dh, cellnum)
    reinit!(c.cellvalues, c.coords)
    _copydofs!(c.ae, anew, c.dofs)     # ae_new = a_new[dofs]
    _copydofs!(c.ae_old, aold, c.dofs) # ae_old = a_old[dofs]
    fill!(c.Ke, 0)
    fill!(c.re, 0)
    return nothing  # Ferrite's reinit! doesn't return 
end

# End of required functions for a custom CellBuffer

"""
    create_threaded_CellBuffers(c::CellBuffer; nthreads=Threads.nthreads())
    create_threaded_CellBuffers(cs::Tuple; nthreads=Threads.nthreads())
    create_threaded_CellBuffers(cs::Dict{String}; nthreads=Threads.nthreads())

Convenience function for creating a vector with cell buffers for each thread. 
The standard workflow is to first call `CellBuffer` with the 
dof handler. For `DofHandler` this will give a `CellBuffer`, 
and for `MixedDofHandler` this gives a tuple of `CellBuffer`s. 
In both cases, the output can be given to `create_threaded_CellBuffers`
to produce the appropriate result required by the threaded versions
of [`doassemble!`](@ref).
Similarily for a grid with mixed materials created by a dictionary of materials, 
just pass the created CellBuffer to this function, and the output is what 
`doassemble!` expects. 
"""
create_threaded_CellBuffers(c; nthreads=Threads.nthreads()) = [deepcopy(c) for _ in 1:nthreads]
function create_threaded_CellBuffers(cs::Tuple; nthreads=Threads.nthreads())
    return map(c->create_threaded_CellBuffers(c;nthreads=nthreads), cs)
end
function create_threaded_CellBuffers(c::Dict{String}; nthreads=Threads.nthreads())
    return Dict(key=>create_threaded_CellBuffers(val; nthreads) for (key,val) in c)
end