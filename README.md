# FerriteAssembly

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://KnutAM.github.io/FerriteAssembly.jl/dev)
[![Build Status](https://github.com/KnutAM/FerriteAssembly.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/KnutAM/FerriteAssembly.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/KnutAM/FerriteAssembly.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/KnutAM/FerriteAssembly.jl)

General assembly system for [Ferrite.jl](https://github.com/Ferrite-FEM/Ferrite.jl/).
1. Define your material struct
2. Overload either `element_routine!` or `element_residual!` for your material
3. Assemble with `doassemble!`

## Selected features
* Threaded assembly
* Efficient autodiff for tangent stiffness
* Multiple materials and fields on subdomains
* State variables

## Using with `Ferrite.jl`'s master branch
`Ferrite#master`: [![Ferrite#master](https://github.com/KnutAM/FerriteAssembly.jl/actions/workflows/FerriteMasterCI.yml/badge.svg?branch=main)](https://github.com/KnutAM/FerriteAssembly.jl/actions/workflows/FerriteMasterCI.yml?query=branch%3Amain)
