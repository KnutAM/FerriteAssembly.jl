```@meta
CurrentModule = FerriteAssembly
```
# Hacking API
The functions described on this page are somewhere inbetween regular API and internal API,
and are intended for making custom solutions, but with a higher risk of breaking changes.
Breaking changes to these interfaces should still be followed by a version bump.

## TaskLocals API (task local storage)
During multithreaded assembly, each task needs its own storage with values that it can change.
This can be both cache variables (whose values don't matter after the task completes) or 
other values that are part of the result from the assembly procedure. To do this in a structured
way, the `TaskLocals` type and associated interface is used for all these cases in the package,
but the user should never "see" this type directly (but it can be embedded in types seen by the user,
such as `ThreadedDomainBuffer`). Specifically, a so-called `scatter-gather` approach is emulated, 
even though the memory is shared, allowing this to be simplified. 
```@docs
FerriteAssembly.TaskLocals
FerriteAssembly.create_local
FerriteAssembly.scatter!
FerriteAssembly.gather!
FerriteAssembly.get_local
FerriteAssembly.get_locals
FerriteAssembly.get_base
```

## Assembler interface
Different types of assemblers can be created in addition to those already defined by the package.
The interface for creating an assembler is that the assembler must support the `TaskLocals` API, 
as well as the method, 
`FerriteAssembly.work_single_cell!(assembler, buffer::AbstractCellBuffer)`,
The following methods for builtin assemblers can be used as examples
```@docs
FerriteAssembly.work_single_cell!
```
In addition, the following functions are used to determine properties/requirements for the assembler. 
- `can_thread(assembler)::Bool`: Is multithreading supported, if `false`, sequential assembly will run even if the a `ThreadedDomainBuffer` is used. If not defined, defaults to `false`. 
- `skip_this_domain(assembler, name::String)::Bool`: When running multiple domains, should the domain with key `name` be skipped? Defaults to `false`. 

## Custom scaling
In general, a scaling must, in addition to the `TaskLocals` API, 
support the following functions
```@docs
FerriteAssembly.update_scaling!
FerriteAssembly.reset_scaling!
```

## Threading model (internal)
```@docs
FerriteAssembly.TaskChunks
```