```@meta
CurrentModule = FerriteAssembly
```

# FerriteAssembly
The goal of [FerriteAssembly](https://github.com/KnutAM/FerriteAssembly.jl) 
is make it convenient to iterate over a grid and perform some task for each item.
The primary use is to iterate over all cells and assemble into global matrices and 
vectors. However, it makes many other tasks easy as well, for example iterating over facets for adding boundary conditions, or integrating a function over a given domain. 

The design is built around two main types of objects
1. **Workers**: What to do for a given item
2. **Domains**: What items to work on

The package user is responsible for writing the code for the actual work (e.g., writing the element routine). A given combination of a worker and a type of domain requires the user to overload a specific function. For example, an assembler (worker) and a cell domain requires overloading [`element_routine!`](@ref) or [`element_residual!`](@ref).

## Documentation structure
The documentation has two main parts

| 1) Learning by doing: *runnable examples*                | 2) Reference: *explanations and docstrings*              | 
| :------------------------------------------------------- | :------------------------------------------------------- |
| Tutorials: Complete examples                             | Domains: *What items to work on*                         |
| Short guides for specific tasks                          | Workers: *What to do for an item*                        | 
|                                                          | Convenience: *Premade solutions to make life easy*       |
|                                                          | Customizations: *When the standard isn't good enough*    |
|                                                          | Internals: *Developer documentation*                     |
