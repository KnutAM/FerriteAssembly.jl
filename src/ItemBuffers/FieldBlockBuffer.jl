"""
    FieldBlockBuffer(buffer::AbstractCellBuffer, n::Int)

View of `buffer` where the local dofs are limited to the first `n` dofs, i.e. to the
fields including their attached algebraic variables, but excluding the standalone algebraic
variables. Created by [`field_block`](@ref FerriteAssembly.field_block).
`get_ae`, `get_aeold`, `get_re`, `get_Ke`, and `celldofs` give views of the first `n` dofs,
and all other accessors are forwarded to `buffer`.
"""
struct FieldBlockBuffer{B <: AbstractCellBuffer} <: AbstractCellBufferWrapper
    buffer::B
    n::Int
end
@inline get_wrapped(b::FieldBlockBuffer) = b.buffer

@inline get_ae(b::FieldBlockBuffer) = view(get_ae(b.buffer), 1:b.n)
@inline get_aeold(b::FieldBlockBuffer) = view(get_aeold(b.buffer), 1:b.n)
@inline get_re(b::FieldBlockBuffer) = view(get_re(b.buffer), 1:b.n)
@inline get_Ke(b::FieldBlockBuffer) = view(get_Ke(b.buffer), 1:b.n, 1:b.n)
@inline Ferrite.celldofs(b::FieldBlockBuffer) = view(celldofs(b.buffer), 1:b.n)

"""
    FerriteAssembly.field_block(Ke, re, ae, buffer::AbstractCellBuffer)

**Experimental**. Get the block of the local degrees of freedom belonging to the fields (including
the degrees of freedom of algebraic variables attached to fields), excluding the standalone algebraic
variables (see the `algebraic_variables` keyword of [`DomainSpec`](@ref)). This makes it possible to
call an existing element routine for the fields from the element routine of a domain that also has
standalone algebraic variables, e.g.
```julia
function FerriteAssembly.element_routine!(Ke, re, state, ae, m::WithMultiplier, cellvalues, buffer)
    Ke_f, re_f, ae_f, buffer_f = FerriteAssembly.field_block(Ke, re, ae, buffer)
    FerriteAssembly.element_routine!(Ke_f, re_f, state, ae_f, m.inner, cellvalues, buffer_f)
    # Add the contributions for the standalone variables, using `dof_range(buffer, name)`
end
```
Returns views `Ke_f`, `re_f`, `ae_f` of the leading block of `Ke`, `re`, and `ae`, and a view
of the buffer, `buffer_f`, for which `get_ae`, `get_aeold`, `get_re`, `get_Ke`, and `celldofs` are limited to
the leading block, while all other accessors give the same as for `buffer`.
`Ke` can be `nothing`, e.g. when called from `element_residual!`, and `Ke_f` is then `nothing`.
The state variables and the cache of the inner routine are those of the outer material, so these should
be forwarded if the inner routine requires them, e.g.
```julia
FerriteAssembly.create_cell_state(m::WithMultiplier, args...) = FerriteAssembly.create_cell_state(m.inner, args...)
FerriteAssembly.allocate_cell_cache(m::WithMultiplier, cv) = FerriteAssembly.allocate_cell_cache(m.inner, cv)
```
For a domain without standalone algebraic variables, the full arrays are viewed.
"""
@inline function field_block(Ke, re::AbstractVector, ae::AbstractVector, buffer::AbstractCellBuffer)
    n = nfield_dofs(get_layout(buffer), length(get_ae(buffer)))
    return _block_view(Ke, n), view(re, 1:n), view(ae, 1:n), FieldBlockBuffer(buffer, n)
end
@inline _block_view(::Nothing, ::Int) = nothing
@inline _block_view(Ke::AbstractMatrix, n::Int) = view(Ke, 1:n, 1:n)

# Number of local dofs of the fields, including the dofs of attached variables
@inline nfield_dofs(::Nothing, ndofs::Int) = ndofs
@inline nfield_dofs(layout::AlgebraicLayout, ::Int) = mapreduce(last, max, layout.field_ranges)
