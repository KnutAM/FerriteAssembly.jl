"""
    DofBuffer(sdh::SubDofHandler)
    DofBuffer(domain::DomainSpec)

The part of an item buffer that follows from the local degrees of freedom (dofs) of an item
(cell or facet): The global dof numbers, `dofs`, the local dof range for each field, `dofrange`,
the local dof values, `ae` and `ae_old`, the local residual, `re`, and the local stiffness, `Ke`.

The `DofBuffer` constructor and [`set_dofs!`](@ref FerriteAssembly.set_dofs!) are the only
functions that define which dofs an item has and how these are ordered locally. The item buffers
and the creation of state variables get this information from a `DofBuffer`.
"""
struct DofBuffer{T, DR}
    ae_old::Vector{T}   # Old local dof values
    ae::Vector{T}       # Current local dof values
    re::Vector{T}       # Residual/force vector
    Ke::Matrix{T}       # Stiffness matrix
    dofs::Vector{Int}   # Global dof numbers
    dofrange::DR        # Local dof range for each field (NamedTuple)
end
function DofBuffer(sdh::SubDofHandler)
    n = ndofs_per_cell(sdh)
    return DofBuffer(zeros(n), zeros(n), zeros(n), zeros(n, n), zeros(Int, n), create_dofrange(sdh))
end

# TaskLocals interface
create_local(d::DofBuffer) = DofBuffer(map(deepcopy, (d.ae_old, d.ae, d.re, d.Ke, d.dofs, d.dofrange))...)

"""
    set_dofs!(d::DofBuffer, dh::DofHandler, cellnum::Int)

Set the global dof numbers in `d` to those for the item belonging to cell number `cellnum`.
"""
@inline set_dofs!(d::DofBuffer, dh, cellnum::Int) = celldofs!(d.dofs, dh, cellnum)

"""
    reinit_buffer!(d::DofBuffer, sim::Simulation, cellnum::Int)

Reinitialize `d` for the item belonging to cell number `cellnum`: Set the global dof numbers with
[`set_dofs!`](@ref FerriteAssembly.set_dofs!), copy the local dof values from the global vectors
`sim.a` (current) and `sim.aold`, and zero the local residual and stiffness. If a global vector is
`nothing`, the corresponding local values are set to `NaN`.
"""
@inline function reinit_buffer!(d::DofBuffer, sim, cellnum::Int)
    set_dofs!(d, get_dofhandler(sim), cellnum)
    _copydofs!(d.ae,     sim.a,    d.dofs) # ae_new .= a_new[dofs]
    _copydofs!(d.ae_old, sim.aold, d.dofs) # ae_old .= a_old[dofs]
    fill!(d.Ke, 0)
    fill!(d.re, 0)
    return nothing
end

"""
    reinit_values!(values, coords, grid, cellnum::Int)
    reinit_values!(values, coords, grid, facet::FacetIndex)

Update the cell coordinates, `coords`, and reinitialize `values` (e.g. `CellValues`
or `FacetValues`, or a `NamedTuple` of such) for the given cell or facet in `grid`.
"""
@inline function reinit_values!(values, coords, grid, cellnum::Int)
    getcoordinates!(coords, grid, cellnum)
    reinit!(values, getcells(grid, cellnum), coords)
    return nothing
end
@inline function reinit_values!(values, coords, grid, (cellnum, facetnr)::FacetIndex)
    getcoordinates!(coords, grid, cellnum)
    reinit!(values, getcells(grid, cellnum), coords, facetnr)
    return nothing
end
