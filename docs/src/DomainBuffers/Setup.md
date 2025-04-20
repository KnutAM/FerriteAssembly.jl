```@meta
CurrentModule = FerriteAssembly
```

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
