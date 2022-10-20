# FerriteAssembly

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://KnutAM.github.io/FerriteAssembly.jl/dev)
[![Build Status](https://github.com/KnutAM/FerriteAssembly.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/KnutAM/FerriteAssembly.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/KnutAM/FerriteAssembly.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/KnutAM/FerriteAssembly.jl)

General assembly system for [Ferrite.jl](https://github.com/Ferrite-FEM/Ferrite.jl/).
Exports `doassemble!`, `CellBuffer`, `create_threaded_CellBuffers`, and `create_threaded_assemblers`.
Requires users to define `element_routine!` or `element_residual!`

## Features
* Less boilerplate code
* Threaded assembly supported
* Element stiffness can be calculated with autodiff
* User-defined scaling of unknowns and residuals during assembly supported

