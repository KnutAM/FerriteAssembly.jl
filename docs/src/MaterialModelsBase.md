# MaterialModelsBase
For constitutive mechanical behavior following the 
[`MaterialModelsBase.jl`](https://github.com/KnutAM/MaterialModelsBase.jl) interface,
[`element_routine!`](@ref FerriteAssembly.element_routine!), 
[`element_residual!`](@ref FerriteAssembly.element_residual!), 
[`create_cell_state`](@ref FerriteAssembly.create_cell_state), 
and [`allocate_cell_cache`](@ref FerriteAssembly.allocate_cell_cache)
are implemented in `FerriteAssembly.jl`.

```@docs
FerriteAssembly.element_routine!(Ke, re, state_new::Vector{<:MMB.AbstractMaterialState}, ae, material::MMB.AbstractMaterial, cellvalues::CellVectorValues, buffer)
FerriteAssembly.element_residual!(re, state_new::Vector{<:MMB.AbstractMaterialState}, ae, material::MMB.AbstractMaterial, cellvalues::CellVectorValues, buffer)
FerriteAssembly.create_cell_state(m::MMB.AbstractMaterial, cv::CellVectorValues, args...)
FerriteAssembly.allocate_cell_cache(m::MMB.AbstractMaterial, ::Any)
```
