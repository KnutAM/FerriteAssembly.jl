# Internal API
Note that the internal API may change without being considered a breaking change!

## Item buffers
An item buffer holds what is required to work one item (cell or facet). The parts are kept in
one place each, such that an extension only needs to change that place:

* The local degrees of freedom (which global dofs an item has, their local order, and the
  arrays with this size) are defined by the `DofBuffer`, which is used by `CellBuffer`,
  `FacetBuffer`, and by the creation of state variables.
* A cell buffer wrapping another cell buffer subtypes `AbstractCellBufferWrapper`, which
  forwards all accessor functions to the wrapped buffer.
* The setup functions take the [`DomainSpec`](@ref) and get the information they need from it.

```@docs
FerriteAssembly.DofBuffer
FerriteAssembly.get_dofbuffer
FerriteAssembly.set_dofs!
FerriteAssembly.CellBuffer
FerriteAssembly.AbstractCellBufferWrapper
FerriteAssembly.AutoDiffCellBuffer
FerriteAssembly.FacetBuffer
FerriteAssembly.CoupledCellBuffer
FerriteAssembly.get_itembuffer
FerriteAssembly.reinit_buffer!
FerriteAssembly.fast_getindex
FerriteAssembly._copydofs!
```

## Algebraic layout
```@docs
FerriteAssembly.AlgebraicLayout
FerriteAssembly.compute_fieldwise_ranges
FerriteAssembly.setup_layout
FerriteAssembly.create_dofrange
FerriteAssembly.fill_celldofs!
FerriteAssembly.update_algebraic_dofs!
```

## States
```@docs
FerriteAssembly.create_states
FerriteAssembly._create_cell_state
```

## DomainBuffers
```@docs
FerriteAssembly._internal_setindex!
FerriteAssembly._internal_deletekeys!
```

## Coupling
```@docs
FerriteAssembly.CoupledSimulation
FerriteAssembly.get_domainbuffer
FerriteAssembly._uncouple
```

## Work execution
```@docs
FerriteAssembly.skip_this_domain
FerriteAssembly.can_thread
FerriteAssembly.work_single_cell!
FerriteAssembly.work_single_facet!
FerriteAssembly.autogenerate_cellvalues
FerriteAssembly.autogenerate_facetvalues
```

## Threading model
```@docs
FerriteAssembly.TaskChunks
FerriteAssembly.drain_tasks
FerriteAssembly.spawn_tasks_safe_wait
```

## Utilities
```@docs
FerriteAssembly.asset_url
```
