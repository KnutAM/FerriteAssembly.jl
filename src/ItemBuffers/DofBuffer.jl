"""
    DofBuffer(sdh::SubDofHandler, layout = nothing)
    DofBuffer(domain::DomainSpec)

The part of an item buffer that follows from the local degrees of freedom (dofs) of an item
(cell or facet): The global dof numbers, `dofs`, the local dof range for each field, `dofrange`,
the local dof values, `ae` and `ae_old`, the local residual, `re`, and the local stiffness, `Ke`.
If the domain uses algebraic variables, `layout` is the [`AlgebraicLayout`](@ref FerriteAssembly.AlgebraicLayout),
which describes how the cell dofs and the dofs of the algebraic variables are ordered locally.

The `DofBuffer` constructor and [`set_dofs!`](@ref FerriteAssembly.set_dofs!) are the only
functions that define which dofs an item has and how these are ordered locally. The item buffers
and the creation of state variables get this information from a `DofBuffer`.
"""
struct DofBuffer{T, DR, L}
    ae_old::Vector{T}   # Old local dof values
    ae::Vector{T}       # Current local dof values
    re::Vector{T}       # Residual/force vector
    Ke::Matrix{T}       # Stiffness matrix
    dofs::Vector{Int}   # Global dof numbers
    scratch::Vector{Int}# Cell dofs, used if the local dofs of a field are not given by `celldofs` directly
    dofrange::DR        # Local dof range for each field and algebraic variable (NamedTuple)
    layout::L           # AlgebraicLayout, or nothing if no algebraic variables are used
end
function DofBuffer(sdh::SubDofHandler, layout = nothing)
    n = ndofs_per_item(sdh, layout)
    dofs = zeros(Int, n)
    update_algebraic_dofs!(dofs, layout, sdh.dh)
    scratch = zeros(Int, scratch_length(layout))
    return DofBuffer(zeros(n), zeros(n), zeros(n), zeros(n, n), dofs, scratch, create_dofrange(sdh, layout), layout)
end

# TaskLocals interface
create_local(d::DofBuffer) = DofBuffer(map(deepcopy, (d.ae_old, d.ae, d.re, d.Ke, d.dofs, d.scratch, d.dofrange))..., d.layout)

"""
    set_dofs!(d::DofBuffer, dh::DofHandler, cellnum::Int)

Set the global dof numbers in `d` to those for the item belonging to cell number `cellnum`.
The dofs of algebraic variables do not depend on the item, and are not changed,
see [`update_algebraic_dofs!`](@ref FerriteAssembly.update_algebraic_dofs!).
"""
@inline set_dofs!(d::DofBuffer, dh, cellnum::Int) = fill_celldofs!(d.dofs, d.scratch, d.layout, dh, cellnum)

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
