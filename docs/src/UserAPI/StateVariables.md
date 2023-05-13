# State variables
The state variable for a given cell is determined by the material type, via 
overloading the [`create_cell_state`](@ref FerriteAssembly.create_cell_state)
function. To update old states to the new states, use [`update_states!`](@ref).
```@docs
FerriteAssembly.create_cell_state(args...)
update_states!
```