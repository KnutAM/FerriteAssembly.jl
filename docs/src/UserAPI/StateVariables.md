## State variables
The initial state variables may vary depending on the position in the grid.
Furthermore, the datastructure depends on the type of dof handler, so
a convenience function exists that creates the correct variable structure. 
```@docs
FerriteAssembly.create_cell_state
update_states!
FerriteAssembly.create_states
```