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
and `cache` to the `CellBuffer` as well. Note that this constructor is normally 
not used, and is instead called from [`setup_cellbuffer`](@ref). 
"""
function CellBuffer(numdofs::Int, numnodes::Int, ::Val{dim}, cellvalues, material, cell_load=nothing, cache=nothing) where dim
    return CellBuffer(
        zeros(numdofs), zeros(numdofs), zeros(numdofs), zeros(numdofs,numdofs), 
        zeros(Int, numdofs), zeros(Vec{dim}, numnodes), 
        cellvalues, material, cell_load, cache
        )
end


"""
    setup_cellbuffer(
        dh::DofHandler, cellvalues, material, 
        cell_load=nothing, cache=nothing)

Creates a single `CellBuffer` for use with the standard `DofHandler` and a single material.
"""
function setup_cellbuffer(dh::DofHandler{dim}, args...) where dim
    return CellBuffer(ndofs_per_cell(dh), Ferrite.nnodes_per_cell(dh), Val{dim}(), args...)
end
"""
    setup_cellbuffer(
        dh::MixedDofHandler, cellvalues, material, 
        cell_load=nothing, cache=nothing)
    
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
function setup_cellbuffer(dh::MixedDofHandler, cellvalues, materials, 
        cell_load=nothing, caches=nothing)
    numfh = length(dh.fieldhandlers)
    cellvalues_ = _maketuple(cellvalues, numfh)
    materials_ = _maketuple(materials, numfh)
    cell_load_ = _maketuple(cell_load, numfh)
    caches_ = _maketuple(caches, numfh)
    return ntuple(i->setup_cellbuffer(dh, dh.fieldhandlers[i], cellvalues_[i], materials_[i], cell_load_[i], caches_[i]), numfh)
end

function setup_cellbuffer(dh::MixedDofHandler{dim}, fh::FieldHandler, args...) where dim
    return CellBuffer(ndofs_per_cell(dh, fh), Ferrite.nnodes_per_cell(dh, fh), Val{dim}(), args...)
end

"""
    setup_cellbuffer(dh::AbstractDofHandler, cv, materials::Dict, 
        cell_load=nothing, cache=nothing)

Return a `Dict{String}` for each `material` in `materials`.
If any of `cv`, `cell_load`, or `cache` is not a `Dict`, the same value 
is used for each `CellBuffer`. If `Dict`s are used, the keys must match 
those in `materials`. The return type depends on the `dh`:

* `dh::DofHandler`: `Dict{String,<:CellBuffer}`
* `dh::MixedDofHandler`: `Dict{String,NTuple{N,CellBuffer}`
"""
function setup_cellbuffer(dh::DofHandler, cv, mtrls::Dict{String}, args...; kwargs...)
    return _setup_dict_cellbuffer(dh, cv, mtrls, args...; kwargs...)
end
function setup_cellbuffer(dh::MixedDofHandler, cv, mtrls::Dict{String}, args...; kwargs...)
    return _setup_dict_cellbuffer(dh, cv, mtrls, args...; kwargs...)
end
# Need special treatment due to type ambiguity
function _setup_dict_cellbuffer(dh::Union{DofHandler,MixedDofHandler}, cv, 
        materials::Dict{String}, cell_load=nothing, caches=nothing)
    setkeys = keys(materials)
    cv_ = _makedict(cv, setkeys)
    cell_load_ = _makedict(cell_load, setkeys)
    caches_ = _makedict(caches, setkeys)
    return Dict(k=>
        setup_cellbuffer(dh, cv_[k], materials[k], cell_load_[k], caches_[k])
        for k in setkeys)
end

# Required functions for a custom CellBuffer (only required internally)
@inline get_Ke(c::CellBuffer) = c.Ke
@inline get_re(c::CellBuffer) = c.re
@inline get_ae(c::CellBuffer) = c.ae
@inline get_material(c::CellBuffer) = c.material
@inline get_cellvalues(c::CellBuffer) = c.cellvalues
@inline Ferrite.celldofs(c::CellBuffer) = c.dofs

# Convenience access functions to CellBuffer for use inside element routines
"""
    Ferrite.getcoordinates(c::CellBuffer)

`Ferrite.jl`'s `getcoordinates` function is overloaded on the `CellBuffer` to return 
the current cell's nodal coordinates. 
"""
@inline Ferrite.getcoordinates(c::CellBuffer) = c.coords

"""
    FerriteAssembly.get_aeold(c::CellBuffer)

Get the old element dof-values for the current cell
"""
@inline get_aeold(c::CellBuffer) = c.ae_old

"""
    FerriteAssembly.get_load(c::CellBuffer)

Get the user specified body load given to `CellBuffer`
"""
@inline get_load(c::CellBuffer) = c.cell_load

"""
    FerriteAssembly.get_cache(c::CellBuffer)

Get the user-specified `cache` given to the `CellBuffer`
"""
@inline get_cache(c::CellBuffer) = c.cache

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