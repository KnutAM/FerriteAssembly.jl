# Integrators
In addition to assembling system matrices and vectors, much of the internal code can be reused 
to create quite efficient integration of values, given a solution vector (and potentially state variables).
While this feature can be used to calculate per-integration point quantities, using the [`QuadPointEvaluator`](@ref) is recommended. 

## Available integrators

### SimpleIntegrator
```@docs
SimpleIntegrator
```

### Integrator
```@docs
Integrator
```

## Functions to overload
Using the `Integrator` requires the following functions 
to be overloaded, depending on the type of domain
```@docs
FerriteAssembly.integrate_cell!
FerriteAssembly.integrate_facet!
```
