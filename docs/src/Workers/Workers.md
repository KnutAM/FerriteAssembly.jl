# Workers
There are currently two categories of workers implemented:
- [Assemblers](@ref): Calculate system matrices and vectors 
- [Integrators](@ref): Integrate a function over the domain

## Perform the work
Having setup a [DomainBuffer](@ref) and a worker, work is performed by calling the function `work!`
```@docs
work!
```

## How a worker works
`work_single_cell!`
`work_single_face!`

### Support threading
`TaskLocals` interface + `can_thread` defined