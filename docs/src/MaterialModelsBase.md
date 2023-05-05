# MaterialModelsBase
If you have implemented your constitutive mechanical behavior following the 
[`MaterialModelsBase.jl`](https://github.com/KnutAM/MaterialModelsBase.jl)'s API,
[`FerriteAssembly.element_routine!`](@ref) and [`FerriteAssembly.create_cell_state`](@ref) 
have already been implemented for this case. 

```@docs
FerriteAssembly.element_routine!(Ke, re, state::Vector{<:FerriteAssembly.MMB.AbstractMaterialState}, ae, material::FerriteAssembly.MMB.AbstractMaterial, cellvalues::CellVectorValues, dh_fh, Δt, cb)
```