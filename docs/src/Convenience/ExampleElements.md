# Example elements
The package includes a set of example elements inside the submodule `FerriteAssembly.ExampleElements`.
In order to use these it is possible to import them explicitly as, e.g.,

`import FerriteAssembly.ExampleElements: StationaryFourier`

### Overview
* [`StationaryFourier`](@ref FerriteAssembly.ExampleElements.StationaryFourier)
* [`TransientFourier`](@ref FerriteAssembly.ExampleElements.TransientFourier)
* [`ElasticPlaneStrain`](@ref FerriteAssembly.ExampleElements.ElasticPlaneStrain)
* [`PoroElasticPlaneStrain`](@ref FerriteAssembly.ExampleElements.PoroElasticPlaneStrain)
* [`WeakForm`](@ref FerriteAssembly.ExampleElements.WeakForm)
* [`J2Plasticity`](@ref FerriteAssembly.ExampleElements.J2Plasticity)

## Available example elements
The following elements are implemented as full elements
```@docs
ExampleElements.StationaryFourier
ExampleElements.TransientFourier
ExampleElements.ElasticPlaneStrain
ExampleElements.PoroElasticPlaneStrain
ExampleElements.WeakForm
```

## Available example mechanical materials
The following physical behaviors are implemented as a 
`MaterialModelsBase.jl` material, which is supported as
a regular element as well.
```@docs 
ExampleElements.J2Plasticity
```