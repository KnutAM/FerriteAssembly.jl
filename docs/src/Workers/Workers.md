# Workers
There are currently three categories of workers implemented:
- [Assemblers](@ref): Calculate system matrices and vectors 
- [Integrators](@ref): Integrate a function over the domain
- [QuadPointEvaluator](@ref): Evaluate a function in each quadrature point

## Perform the work
Having setup a [`DomainBuffer`](@ref Setup-API) and a worker, work is performed by calling the function `work!`
```@docs
work!
```

A worker needs to define what to do for different domains, 
specifically it needs to define 
`work_single_cell!` and `work_single_facet!`. Further details 
about writing custom workers can be found in the 
[Customizations](@ref Assembler-interface) section.