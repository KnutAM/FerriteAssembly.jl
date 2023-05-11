## Residual scaling
There are many options for how to scale the residual in finite element simulations.
This package does not intend to implement many different options, but does give the 
user the option to calculate scaling contributions from each cell, which may be useful.
By defining a `scaling` that is passed to [`setup_assembly`](@ref), it can be updated 
based on the output from each cell. 

One type of scaling, [`ElementResidualScaling`](@ref), is included as described below.
Its code can be used as a template for how to include custom scaling that works on 
the element level.

### ElementResidualScaling
```@docs
ElementResidualScaling
```
