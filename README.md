# FerriteAssembly

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://KnutAM.github.io/FerriteAssembly.jl/dev)
[![Build Status](https://github.com/KnutAM/FerriteAssembly.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/KnutAM/FerriteAssembly.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/KnutAM/FerriteAssembly.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/KnutAM/FerriteAssembly.jl)

# FerriteAssembly
The goal of [FerriteAssembly](https://github.com/KnutAM/FerriteAssembly.jl) 
is to provide a simple, but versatile and efficient, structure for assembling in 
[Ferrite.jl](https://github.com/Ferrite-FEM/Ferrite.jl/).

At its core, `FerriteAssembly.jl` lets you overload some of the following functions
```julia
element_routine!(Ke, re, state, ae, material, cellvalues, cellbuffer)
element_residual!(re, state, ae, material, cellvalues, cellbuffer)
face_routine!(Ke, re, ae, material, facevalues, facebuffer)
face_residual!(re, ae, material, facevalues, facebuffer)
```
define a `domainbuffer` with `setup_domainbuffer`,

and then call `work!(assembler, domainbuffer)` 
with an `assembler` of your choice. 

This makes it easy to share different elements between researchers,
and provides a structured way to solve various `Ferrite.jl` problems. 

## Key features
* **What FerriteAssembly helps with**
  * Assemble into system vectors and/or matrices
  * Integrate functions over the domain
  * Postprocess cell data over the domain
* **Convenience features**
  * `LoadHandler` for adding Neumann boundary conditions and body loads/source terms as easy as Dirichlet boundary conditions with `Ferrite`'s `ConstraintHandler`.
  * `ExampleElements` submodule with various behaviors for quick testing
  * Built-in support for mechanical materials following the `MaterialModelsBase.jl` interface
* **Supports many different cases**
  * Easy switching between *sequential* and *threaded* work.
  * *Multiple domains* with different fields, interpolations, and/or element routines.
  * Efficient *automatic differentiation* for tangent stiffness. 
  * Handling old and new *state variables*

See the [documentation](https://KnutAM.github.io/FerriteAssembly.jl/dev) for more details.

## Using with `Ferrite.jl`'s master branch
This badge shows if `FerriteAssembly#main` is compatible with `Ferrite#master`.
Currently, it will likely remain incompatible until some open issues regarding interfaces for the new Ferrite v0.4 release are solved. 

[![Ferrite#master](https://github.com/KnutAM/FerriteAssembly.jl/actions/workflows/FerriteMasterCI.yml/badge.svg?branch=main)](https://github.com/KnutAM/FerriteAssembly.jl/actions/workflows/FerriteMasterCI.yml?query=branch%3Amain)
