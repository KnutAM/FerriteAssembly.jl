# State variables
The state variable for a given cell is determined by the material type, via 
overloading the [`create_cell_state`](@ref FerriteAssembly.create_cell_state)
function. To update old states to the new states, use [`update_states!`](@ref).

## The state variable datastructure
The state variables are created when calling [`setup_assembly`](@ref).
For a single domain, the output is `Dict{Int}`, where each element is the 
output of `create_cell_state`. For multiple domains, the output is 
`Dict{String,Dict{Int}}` with keys from the `AssemblyDomain`s, and the 
contained `Dict{Int}` contain the output from `create_cell_state`. 
In all cases, the `Dict{Int}`s are indexed by the cell nr. 

## API
```@docs
FerriteAssembly.create_cell_state(args...)
update_states!
```