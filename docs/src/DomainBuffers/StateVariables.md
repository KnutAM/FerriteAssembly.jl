# State variables
The state variable for a given cell is determined by the material type, via 
overloading the [`create_cell_state`](@ref FerriteAssembly.create_cell_state)
function. To update old states to the new (just-converged) states, use
[`update_states!`](@ref update_states!(::FerriteAssembly.DomainBuffers)), which by default
copies the values so that both old and new states hold the converged values
directly afterwards (safe to read immediately, e.g. for postprocessing).

## Special cases
Two typical cases are considered by default. A type, `CS`, which defines the state
for the entire cell and a cell state that consists of an `AbstractVector{QS}` where the type
`QS` is the type of the state for each quadrature point. In the following both cases are 
described considering the type of the state `CS` or `QS` being denoted `S` as the type of 
the state.

### Non isbits state
If `S` is not a bits type (`isbitstype(S) = false`), then the user must define either 
[`FerriteAssembly.copy_state`](@ref) or a [`FerriteAssembly.copy_state!`](@ref) for 
the default [`update_states!`](@ref update_states!(::FerriteAssembly.DomainBuffers)) to work.

Alternatively, `update_states!` can be called with the keyword argument `mode = :flip` to just
flip the references. However, this comes with the important caveat that one should **never** 
read values from the current state, see details in 
[`update_states!`](@ref update_states!(::FerriteAssembly.DomainBuffers)). This is useful for 
cases when the state consists of large data structures such as with ``\mathrm{FE}^2`` simulations.

### Reverting the states
If the current (new) state is used in the element routines (e.g. as initial guess),
the function [`revert_states!`](@ref revert_states!(::FerriteAssembly.DomainBuffers))
can be used to update the states such that `states = old_states` before retrying to find 
the solution after a failed time step

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
FerriteAssembly.copy_state
FerriteAssembly.copy_state!
FerriteAssembly.remove_dual
```