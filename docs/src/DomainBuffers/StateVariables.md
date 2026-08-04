# State variables
The state variable for a given cell is determined by the material type, via 
overloading the [`create_cell_state`](@ref FerriteAssembly.create_cell_state)
function. To update old states to the new states, use [`update_states!`](@ref update_states!(::FerriteAssembly.DomainBuffers)).
To instead reset the new states back to the old (converged) states, e.g. when retrying
a non-converged increment, use [`set_new_to_old_states!`](@ref set_new_to_old_states!(::FerriteAssembly.DomainBuffers)).
If [`create_cell_state`](@ref FerriteAssembly.create_cell_state) returns a *mutable*
`AbstractArray` (which must keep the same axes between calls, otherwise an `ArgumentError`
is thrown), each element is copied individually; otherwise (including an immutable
`AbstractArray`, e.g. built from `NTuple`s or `StaticArrays`) the whole cell state is
copied. `isbits` values are copied by identity (no allocation); any other value must have a
[`FerriteAssembly.copy_state`](@ref) method for its type — otherwise `set_new_to_old_states!`
throws a `MethodError`.

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
FerriteAssembly.copy_state
FerriteAssembly.remove_dual
```