# Internal API
Note that the internal API may change without being considered a breaking change!

## Item buffers
```@docs
FerriteAssembly.CellBuffer
FerriteAssembly.AutoDiffCellBuffer
FerriteAssembly.FacetBuffer
FerriteAssembly.CoupledCellBuffer
FerriteAssembly.get_itembuffer
FerriteAssembly.reinit_buffer!
FerriteAssembly.fast_getindex
FerriteAssembly._copydofs!
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
