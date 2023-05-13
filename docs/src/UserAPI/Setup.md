```@meta
CurrentModule = FerriteAssembly
```

# Setup
* To setup the assembly, call [`setup_assembly`](@ref)
* For multiple domains, [`AssemblyDomain`](@ref)s must be created first.
* Details about the [`buffer`](@ref AbstractDomainBuffer) output. 
* State variables discussed [here](@ref State-variables)

## Setup API
```@docs
setup_assembly
AssemblyDomain
```

## AbstractDomainBuffer
The `buffer` output from `setup_assembly` will be of type `DomainBuffer`, `ThreadedDomainBuffer`, 
or a `Dict{String}` with eltype of one of the former. A few convenience functions are defined 
for these output types
```@docs
FerriteAssembly.get_material(buffer::DomainBuffer)
```
