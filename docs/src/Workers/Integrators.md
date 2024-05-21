# Integrators
In addition to assembling system matrices and vectors, much of the internal code can be reused 
to create quite efficient integration of values, given a solution vector (and potentially state variables).
This feature can also be used to calculate per-cell (or per-integration point) quantities, simply by 
defining an integrator containing a vector that is indexed by the global cell id. 

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
