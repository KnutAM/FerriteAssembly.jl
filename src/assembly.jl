""" 
    doassemble!(
        K::AbstractMatrix, r::AbstractVector, a::AbstractVector, 
        aold::AbstractVector, s::AbstractVector, dh::DofHandler, 
        cellvalues, material, Δt::Number, cache::CellCache, 
        cellset=nothing
        )

Assemble all cells by using the regular DofHandler `dh`. 
`K`, and `r` are the global stiffness matrix and residual 
vector to be calculated. 
`a`, and `aold`are the current and old unknowns. 
`s` is a vector of state variables, where each
element contains state variables for each cell. 
`cellvalues` are passed on to the element_routine after reinitialization. 
For multiple fields, it can be a Tuple or NamedTuple. 
`material` is a user defined type passed onto the element, 
and so is the time increment `Δt`. 
`cache` is described in [`CellCache`](@ref).
"""
function doassemble!(
    K::AbstractMatrix, r::AbstractVector, a::AbstractVector, 
    aold::AbstractVector, s::AbstractVector, dh::DofHandler, 
    cellvalues, material, Δt::Number, cache::CellCache, 
    cellset=nothing
    )
    assembler = start_assemble(K, r)
    for (i, cell) in enumerate(CellIterator(dh, cellset))
        assemble_cell!(
            assembler, cell, cellvalues, material, 
            s[i], a, aold, dh, Δt, cache)
    end
end

""" doassemble!(
    K::AbstractMatrix, r::AbstractVector, anew::AbstractVector, 
    aold::AbstractVector, state::Tuple, dh::MixedDofHandler, 
    cellvalues::Tuple, materials::Tuple, Δt::Number, caches::Tuple,
    fullcellset=nothing
    )

Outer assembling loop using the MixedDofHandler `dh`
"""
function doassemble!(
    K::AbstractMatrix, r::AbstractVector, anew::AbstractVector, 
    aold::AbstractVector, state::Tuple, dh::MixedDofHandler, 
    cellvalues::Tuple, materials::Tuple, Δt::Number, caches::Tuple,
    fullcellset=nothing)
    assembler = start_assemble(K, r)
    for (fh, cv, m, s, c) in zip(dh.fieldhandlers, cellvalues, materials, state, caches)
        doassemble!(assembler, fh, cv, m, s, dh, anew, aold, Δt, c, fullcellset)
    end
end

""" 
    doassemble!(
        assembler, fh::FieldHandler, cellvalues, material, 
        s::AbstractVector, dh::MixedDofHandler, 
        anew::AbstractVector, aold::AbstractVector, 
        Δt::Number, cache::CellCache, fullcellset
        )

Inner assembling loop using the MixedDofHandler `dh` for a specific field `fh`.
Normally not called by the user, only called from the other
`doassemble!(..., ::MixedDofHandler, ...)`
"""
function doassemble!(
    assembler, fh::FieldHandler, cellvalues, material, 
    s::AbstractVector, dh::MixedDofHandler, 
    anew::AbstractVector, aold::AbstractVector, 
    Δt::Number, cache::CellCache, fullcellset
    )
    if isnothing(fullcellset)
        cellset = collect(fh.cellset)
    else
        cellset = collect(union(fh.cellset, fullcellset))
    end
    for (i, cell) in enumerate(CellIterator(dh, cellset))
        assemble_cell!(
            assembler, cell, cellvalues, material, s[i], 
            anew, aold, fh, Δt, cache)
    end
end

"""
    assemble_cell!(
        assembler, cell, cellvalues, material, 
        state, anew, aold, dh_fh, Δt, cellcache)

Assemble the specific cell, provides a function barrier and uniform 
initialization of cellvalues, zero out element stiffness and residual,
and scaling of primary and residual values
"""
function assemble_cell!(
    assembler, cell, cellvalues, material, 
    state, anew, aold, dh_fh, Δt, cellcache::CellCache
    )
    reinit!(cellcache, cell, anew, aold)
    reinit!(cellvalues, cell)
    ae_new, ae_old, re, Ke, materialcache = getcellcachevars(cellcache)

    unscale_primary!.((ae_new, ae_old), (material,), (dh_fh,))
    element_routine!(
        Ke, re, ae_new, ae_old, state, material, 
        cellvalues, dh_fh, Δt, materialcache
        )
    scale_residual!(Ke, re, material, dh_fh)

    assemble!(assembler, celldofs(cell), Ke, re)
end

Ferrite.reinit!(cvs::Union{Tuple,NamedTuple}, cell::CellIterator) = reinit!.(cvs, (cell,))


