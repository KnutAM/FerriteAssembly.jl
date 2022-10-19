""" 
    doassemble!(
        K::AbstractMatrix, r::AbstractVector, a::AbstractVector, 
        aold::AbstractVector, s::AbstractVector, dh::DofHandler, 
        cellvalues, material, Δt::Number, cache::CellCache
        )

Assemble all cells by using `dh::DofHandler`.
* `K`, and `r` are the global stiffness matrix and residual 
  vector to be calculated. 
* `a` and `aold` are the current and old unknowns. 
* `s` is a vector of state variables, where each
  element contains state variables for each `cellid`. 
  If `cellset` is given (see below), the number is still the `cellid`,
  and in some cases it might make sense to pass a sparse vector. 
* `cellvalues` are passed on to the element routine 
  after calling `Ferrite.reinit!`. (Multiple cellvalues 
  can be given as a `Tuple` or `NamedTuple`, and all
  elements will be reinitialized)
* The user-defined `material` and the time increment `Δt`
  are passed into the element routine. 
* `cache` is described in [`CellCache`](@ref).
"""
function doassemble!(
    K::AbstractMatrix, r::AbstractVector, a::AbstractVector, 
    aold::AbstractVector, s::AbstractVector, dh::DofHandler, 
    cellvalues, material, Δt::Number, cache::CellCache
    )
    assembler = start_assemble(K, r)
    for cell in CellIterator(dh)
        assemble_cell!(
            assembler, cell, cellvalues, material, 
            s[cellid(cell)], a, aold, dh, Δt, cache)
    end
end

""" 
    doassemble!(
        K::AbstractMatrix, r::AbstractVector, anew::AbstractVector, 
        aold::AbstractVector, state::Tuple, dh::MixedDofHandler, 
        cellvalues::Tuple, materials::Tuple, Δt::Number, caches::Tuple
        )

Assemble all cells by using `dh::MixedDofHandler`. In this case, 
we first loop over each `fieldhandler` in `dh`. Therefore, 
some variables must now be given as tuples, where each element corresponds 
to the values in the specific `fieldhandler`. 
* `K`, and `r` are the global stiffness matrix and residual 
  vector to be calculated. 
* `a` and `aold` are the current and old unknowns. 
* `s` is a tuple of vectors of state variables. Each element, `s[i]`,
  corresponds to the `i`th element in `dh.fieldhandlers`. 
  The `j`th element of `s[i]::AbstractVector`, contains the state variables 
  for `cellid=j`. If the state-type differs for multiple fieldhandlers, it 
  might make sense to use sparse vectors. 
* `cellvalues` is a tuple with one element for each fieldhandler. 
  For each fieldhandler `i`, `cellvalues[i]` is passed to the element 
  routine after calling `Ferrite.reinit!`. (Multiple cellvalues can 
  be given by making `cellvalues[i]` a `Tuple` or `NamedTuple`. 
  All elements will be reinitialized)
* The user-defined `material` and the time increment `Δt`
  are passed into the element routine. Finally,
* `caches` is a tuple with `CellCache`-elements, which are 
  described in [`CellCache`](@ref).
"""
function doassemble!(
    K::AbstractMatrix, r::AbstractVector, a::AbstractVector, 
    aold::AbstractVector, state::Tuple, dh::MixedDofHandler, 
    cellvalues::Tuple, materials::Tuple, Δt::Number, caches::Tuple
    )
    assembler = start_assemble(K, r)
    for (fh, cv, m, s, c) in zip(dh.fieldhandlers, cellvalues, materials, state, caches)
        inner_doassemble!(assembler, fh, cv, m, s, dh, a, aold, Δt, c)
    end
end

""" 
    inner_doassemble!(
        assembler, fh::FieldHandler, cellvalues, material, 
        s::AbstractVector, dh::MixedDofHandler, 
        anew::AbstractVector, aold::AbstractVector, 
        Δt::Number, cache::CellCache
        )

Inner assembling loop using the MixedDofHandler `dh` for a specific field `fh`.
Normally not called by the user, only called from
`doassemble!(K, r, a, ..., dh::MixedDofHandler, ...)`
"""
function inner_doassemble!(
    assembler, fh::FieldHandler, cellvalues, material, 
    s::AbstractVector, dh::MixedDofHandler, 
    anew::AbstractVector, aold::AbstractVector, 
    Δt::Number, cache::CellCache
    )
    for cell in CellIterator(dh, collect(fh.cellset))
        assemble_cell!(
            assembler, cell, cellvalues, material, s[cellid(cell)], 
            anew, aold, fh, Δt, cache)
    end
end

# Make reinit! work for Tuple/NamedTuple of FEValues (cellvalues, facevalues)
Ferrite.reinit!(fevalues::Union{Tuple,NamedTuple}, args...) = foreach(fevalue->reinit!(fevalue, args...), fevalues)

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
    (;ae, ae_old, re, Ke, materialcache) = getcellcachevars(cellcache)

    unscale_primary!.((ae, ae_old), (material,), (dh_fh,))
    element_routine!(
        Ke, re, ae, ae_old, state, material, 
        cellvalues, dh_fh, Δt, materialcache
        )
    scale_residual!(Ke, re, material, dh_fh)

    assemble!(assembler, celldofs(cell), Ke, re)
end