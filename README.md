# FerriteAssembly

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://KnutAM.github.io/FerriteAssembly.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://KnutAM.github.io/FerriteAssembly.jl/dev)
[![Build Status](https://github.com/KnutAM/FerriteAssembly.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/KnutAM/FerriteAssembly.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/KnutAM/FerriteAssembly.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/KnutAM/FerriteAssembly.jl)

General and performant assembly system for [Ferrite.jl](https://github.com/Ferrite-FEM/Ferrite.jl/).
Currently only sequential assembly supported, but this will be extended 
to parallel assembly. Both `DofHandler` and `MixedDofHandler` are supported.

Exports `doassemble!` and requires users to define `element_routine!` or `element_residual!`

## Features
* Less boilerplate code
* Element stiffness can be calculated with autodiff
* Scaling of unknowns and residuals

## Future plans
* Support Neumann boundary conditions
* Parallel assembly
