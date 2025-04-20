# Quadrature point data evaluation
While an [`Integrator`](@ref) can be used to evaluate quadrature point data, the [`QuadPointEvaluator`](@ref) provides a more convenient and efficient way to do so. This evaluation is particularly convenient for use with `Ferrite`'s `L2Projector`. 

## `QuadPointEvaluator`
```@docs
QuadPointEvaluator
```

## Functions to overload
For custom quadrature evaluation, `eval_quadpoints_cell!` can be overloaded.
```@docs
FerriteAssembly.eval_quadpoints_cell!
```
