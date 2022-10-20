struct CellBuffer{dim,T,CV,MT,CT}
    ae_old::Vector{T}           # Old element dof values
    ae::Vector{T}               # New element dof values
    re::Vector{T}               # Residual/force vector 
    Ke::Matrix{T}               # Element stiffness matrix
    dofs::Vector{Int}           # celldofs
    coords::Vector{Vec{dim,T}}  # cellcoords
    cellvalues::CV              # Can also be e.g. tuple of cellvalues
    material::MT                # User material definition
    cache::CT                   # CellBuffer as user pleases
end

"""
    CellBuffer(numdofs::Int, numnodes::Int, ::Val{dim}, cellvalues, material, cache=nothing)

Create a cell cache for an element with `numdofs` degrees of freedom and
`numnodes` nodes with dimension `dim`. Add the given `cellvalues`, `material`, 
and `cache` to the `CellBuffer` as well. 

    CellBuffer(dh::DofHandler, cellvalues, material, cache=nothing)

Use `dh` to get `numdofs`, `numnodes`, and `dim`, before calling the above method definition. 

    CellBuffer(dh::MixedDofHandler, fh::FieldHandler, cellvalues, material, cache=nothing)

Use `dh` and `fh` to get `numdofs`, `numnodes`, and `dim`, 
before calling the first `CellBuffer` method definition. 

    CellBuffer(dh::MixedDofHandler, cellvalues, material, cache=nothing)
    
Return a tuple of `CellBuffer`s for each `FieldHandler` in `dh.fieldhandlers`.
`cellvalues[i]` corresponds to `dh.fieldhandlers[i]`, and so does 
`materials[i]` and `caches[i]`. If only one `CellValues`, `material` (not a tuple), and/or `cache`
is given, the same is used for all `fieldhandlers`. If a tuple of `cellvalues` (or materials/caches) 
should be used for each cell, 
and the same tuple should be used for each fieldhandler, 
then it must be given as a tuple of tuples. 
(Alternatively, it is recommended to use a `NamedTuple`)
"""
function CellBuffer(numdofs::Int, numnodes::Int, ::Val{dim}, cellvalues, material, cache=nothing) where dim
    return CellBuffer(
        zeros(numdofs), zeros(numdofs), zeros(numdofs), zeros(numdofs,numdofs), 
        zeros(Int, numdofs), zeros(Vec{dim}, numnodes), 
        cellvalues, material, cache
        )
end

function CellBuffer(dh::DofHandler{dim}, args...) where dim
    return CellBuffer(ndofs_per_cell(dh), Ferrite.nnodes_per_cell(dh), Val{dim}(), args...)
end

function CellBuffer(dh::MixedDofHandler{dim}, fh::FieldHandler, args...) where dim
    return CellBuffer(ndofs_per_cell(dh, fh), Ferrite.nnodes_per_cell(dh, fh), Val{dim}(), args...)
end

function CellBuffer(dh::MixedDofHandler, cellvalues::Union{Tuple,CellValues}, materials, caches=nothing)
    numfh = length(dh.fieldhandlers)
    cellvalues_ = _maketuple(cellvalues, numfh)
    materials_ = _maketuple(materials, numfh)
    caches_ = _maketuple(caches, numfh)
    return ntuple(i->CellBuffer(dh, dh.fieldhandlers[i], cellvalues_[i], materials_[i], caches_[i]), numfh)
end

"""
    Ferrite.reinit!(c::CellBuffer, dh::AbstractDofHandler, cellnum::Int, anew, aold)

Reinitialize the `c::CellBuffer` for cell number `cellnum`.
The global degree of freedom vectors `anew` (current) and `aold` are used
to update the cell degree of freedom vectors in `c`.
The element stiffness, `c.Ke`, and residual, `c.re`, are also zeroed. 
"""
function Ferrite.reinit!(c::CellBuffer, dh::Ferrite.AbstractDofHandler, cellnum::Int, anew, aold)
    celldofs!(c.dofs, dh, cellnum)
    Ferrite.cellcoords!(c.coords, dh, cellnum)
    reinit!(c.cellvalues, c.coords)
    _copyto!(c.ae, anew, c.dofs)     # ae_new = a_new[dofs]
    _copyto!(c.ae_old, aold, c.dofs) # ae_old = a_old[dofs]
    fill!(c.Ke, 0)
    fill!(c.re, 0)
    return nothing  # Ferrite's reinit! doesn't return 
end

"""
    create_threaded_CellBuffers(c::CellBuffer; nthreads=Threads.nthreads())
    create_threaded_CellBuffers(cs::Tuple; nthreads=Threads.nthreads())

Convenience function for creating cell buffers for each thread. 
The standard workflow is to first call `CellBuffer` with the 
dof handler. For `DofHandler` this will give a `CellBuffer`, 
and for `MixedDofHandler` this gives a tuple of `CellBuffer`s. 
In both cases, the output can be given to `create_threaded_CellBuffers`
to produce the appropriate result required by the threaded versions
of [`doassemble!`](@ref).
"""
create_threaded_CellBuffers(c::CellBuffer; nthreads=Threads.nthreads()) = [deepcopy(c) for _ in 1:nthreads]
function create_threaded_CellBuffers(cs::Tuple; nthreads=Threads.nthreads())
    return map(c->create_threaded_CellBuffers(c;nthreads=nthreads), cs)
end