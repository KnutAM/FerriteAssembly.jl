# Internal API
Note that the internal API may change without being considered a breaking change!

```@docs
FerriteAssembly.reinit_buffer!
FerriteAssembly.create_states
FerriteAssembly._copydofs!
FerriteAssembly._create_cell_state
FerriteAssembly.CellBuffer
FerriteAssembly.AutoDiffCellBuffer
FerriteAssembly.FaceBuffer
FerriteAssembly.get_itembuffer
FerriteAssembly.skip_this_domain
FerriteAssembly.can_thread
FerriteAssembly.fast_getindex
FerriteAssembly.work_single_cell!
FerriteAssembly.work_single_face!
FerriteAssembly.autogenerate_cellvalues
FerriteAssembly.autogenerate_facevalues
```

## Threading model
```@docs
FerriteAssembly.TaskChunks
```