# State variables
The state variable for a given cell is determined by the material type, via 
overloading the [`create_cell_state`](@ref FerriteAssembly.create_cell_state)
function. To update old states to the new (just-converged) states, use
[`update_states!`](@ref update_states!(::FerriteAssembly.DomainBuffers)), which by default
copies the values so that both old and new states correctly hold the converged values
directly afterwards (safe to read immediately, e.g. for postprocessing). Pass
`mode = :flip` for the cheaper, allocation-free reference-swap behavior instead — see that
docstring for the gotcha this introduces.

To instead reset the new states back to the old (converged) states, e.g. when retrying
a non-converged increment, use
[`revert_states!`](@ref revert_states!(::FerriteAssembly.DomainBuffers))
(`set_new_to_old_states!` is a deprecated alias for this function).
If [`create_cell_state`](@ref FerriteAssembly.create_cell_state) returns a *mutable*
`AbstractArray` (which must keep the same axes between calls, otherwise an `ArgumentError`
is thrown), each element is copied individually; otherwise (including an immutable
`AbstractArray`, e.g. built from `NTuple`s or `StaticArrays`) the whole cell state is
copied. `isbits` values are copied by identity (no allocation); any other value must have
either a [`FerriteAssembly.copy_state`](@ref) or a [`FerriteAssembly.copy_state!`](@ref)
method for its type — neither has a default method, so a type with neither overloaded
throws a `MethodError`.

[`FerriteAssembly.copy_state!`](@ref) is the allocation-avoiding alternative: for a value
that is itself immutable but wraps a mutable payload (e.g.
`struct MyState; vals::Vector{Float64}; end`), it overwrites the existing destination value
in place (e.g. via `copyto!`) rather than allocating a full replacement, and takes
precedence over [`FerriteAssembly.copy_state`](@ref) when both are applicable. A value's
type needs at most one of the two overloaded, never both.

## The state variable datastructure
The state variables are created when calling [`setup_domainbuffer`](@ref)
or [`setup_domainbuffers`](@ref), and stored inside the buffers. 
The states for a given domain are accessed with 
[`get_state`](@ref FerriteAssembly.get_state(::FerriteAssembly.DomainBuffers, ::String)) 
and 
[`get_old_state`](@ref FerriteAssembly.get_state(::FerriteAssembly.DomainBuffers, ::String)), 
where the state for a particular cell is indexed
by its cell number. (The output from the mentioned functions are `Dict{Int}`)

## API
```@docs
FerriteAssembly.create_cell_state
update_states!
revert_states!
set_new_to_old_states!
FerriteAssembly.copy_state
FerriteAssembly.copy_state!
FerriteAssembly.remove_dual
```