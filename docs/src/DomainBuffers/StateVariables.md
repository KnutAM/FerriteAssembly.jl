# State variables
The state variable for a given cell is determined by the material type, via 
overloading the [`create_cell_state`](@ref FerriteAssembly.create_cell_state)
function. To update old states to the new (just-converged) states, use
[`update_states!`](@ref update_states!(::FerriteAssembly.DomainBuffers)), which by default
copies the values so that both old and new states correctly hold the converged values
directly afterwards (safe to read immediately, e.g. for postprocessing). Pass
`mode = :flip` for the cheaper, allocation-free reference-swap behavior instead, but note
that afterwards the *new* states hold the values from *before* this step until the next
`work!` call — see the gotcha described in that docstring, and in the
[multiple materials](@ref Multiple-materials) and [phase-field fracture](@ref Phase-field-fracture)
tutorials.

To instead reset the new states back to the old (converged) states, e.g. when retrying
a non-converged increment, use [`set_new_to_old_states!`](@ref set_new_to_old_states!(::FerriteAssembly.DomainBuffers)).
The opposite copy direction, updating the old states from the new ones without swapping
references (what `update_states!`'s default `mode = :copy` uses internally), is available as
[`set_old_to_new_states!`](@ref set_old_to_new_states!(::FerriteAssembly.DomainBuffers)).

For both of these copying functions, if [`create_cell_state`](@ref FerriteAssembly.create_cell_state)
returns a *mutable* `AbstractArray` (which must keep the same axes between calls, otherwise an
`ArgumentError` is thrown), each element is copied individually; otherwise (including an
immutable `AbstractArray`, e.g. built from `NTuple`s or `StaticArrays`) the whole cell state
is copied. `isbits` values are copied by identity (no allocation); any other value must have a
[`FerriteAssembly.copy_state`](@ref) method for its type — otherwise a `MethodError` is
thrown.

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
set_new_to_old_states!
set_old_to_new_states!
FerriteAssembly.copy_state
FerriteAssembly.remove_dual
```