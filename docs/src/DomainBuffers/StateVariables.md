# State variables
The state variable for a given cell is determined by the material type, via 
overloading the [`create_cell_state`](@ref FerriteAssembly.create_cell_state)
function. To update old states to the new states, use [`update_states!`](@ref update_states!(::FerriteAssembly.DomainBuffers)).

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
FerriteAssembly.remove_dual
```