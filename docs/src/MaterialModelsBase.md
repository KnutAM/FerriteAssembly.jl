# MaterialModelsBase
If you have implemented your constitutive mechanical behavior following the 
[`MaterialModelsBase.jl`](https://github.com/KnutAM/MaterialModelsBase.jl)'s API,
[`FerriteAssembly.element_routine!`](@ref) and [`FerriteAssembly.create_cell_state`](@ref) 
have already been implemented for this case. 

```@docs
FerriteAssembly.element_routine!(Ke, re, state_new::Vector{<:MMB.AbstractMaterialState}, ae, material::MMB.AbstractMaterial, cellvalues::CellVectorValues, buffer)
```