# FerriteAssembly

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://KnutAM.github.io/FerriteAssembly.jl/dev)
[![Build Status](https://github.com/KnutAM/FerriteAssembly.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/KnutAM/FerriteAssembly.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/KnutAM/FerriteAssembly.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/KnutAM/FerriteAssembly.jl)

# FerriteAssembly
The goal of [FerriteAssembly](https://github.com/KnutAM/FerriteAssembly.jl) 
is to provide a simple, but versatile and efficient, structure for assembling in 
[Ferrite.jl](https://github.com/Ferrite-FEM/Ferrite.jl/).

## Key features
* Easy switching between *sequential* and *threaded* assembly
* *Multiple domains* with different fields, interpolations, and/or element routines.
* Efficient *automatic differentiation* if analytical tangent is not implemented. 
* Support for handling of (old and new) *state variables*
* Easy *integration* of a function over the domain

See the [documentation](https://KnutAM.github.io/FerriteAssembly.jl/dev) for more details.

## Typical workflow
1. Define your custom type and associated element routine (See [Example elements](@ref))
2. Setup Ferrite's `DofHandler` and `CellValues` as usual, and call [`setup_assembly`](@ref)
3. For each assembly, call [`doassemble!`](@ref)

## Using with `Ferrite.jl`'s master branch
This badge shows if `FerriteAssembly#main` is compatible with `Ferrite#master`.
Currently, it will likely remain incompatible until some open issues regarding interfaces for the new Ferrite v0.4 release are solved. 

[![Ferrite#master](https://github.com/KnutAM/FerriteAssembly.jl/actions/workflows/FerriteMasterCI.yml/badge.svg?branch=main)](https://github.com/KnutAM/FerriteAssembly.jl/actions/workflows/FerriteMasterCI.yml?query=branch%3Amain)
