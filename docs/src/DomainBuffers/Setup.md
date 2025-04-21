```@meta
CurrentModule = FerriteAssembly
```
# [Setup](@id DomainBufferSetup)
## Setup API
```@docs
DomainSpec
setup_domainbuffer
setup_domainbuffers
```

## AbstractDomainBuffer
The domain buffer be a `DomainBuffer`, `ThreadedDomainBuffer`, or a `Dict{String}` with eltype 
of one of the former. The following functions are defined for these buffers:
```@docs
FerriteAssembly.get_material(::FerriteAssembly.DomainBuffers, ::String)
FerriteAssembly.get_dofhandler(::FerriteAssembly.DomainBuffers)
FerriteAssembly.get_state(::FerriteAssembly.DomainBuffers, ::String)
FerriteAssembly.get_old_state(::FerriteAssembly.DomainBuffers, ::String)
FerriteAssembly.getset
update_states!(::FerriteAssembly.DomainBuffers)
set_time_increment!(::FerriteAssembly.DomainBuffers, ::Any)
```

## Coupled simulations
The `Simulation` type contains an abstract domain buffer, along with (optionally) 
the global degree of freedom values, which are used to get the local values for each item.
The main purpose is to conveniently collect these when passing into [`work!`](@ref), 
especially in the case of `CoupledSimulations`.

The idea behind the coupled simulation setup is to give access to values from a different simulation
at the item level. For example, when solving two separate problems in parallel, and using staggered
iterations. See the [Phase-field fracture tutorial](@ref Phase-field-fracture) for an example. 
```@docs
Simulation
couple_buffers
CoupledSimulations
```
