```@meta
CurrentModule = FerriteAssembly
```

# FerriteAssembly
The goal of [FerriteAssembly](https://github.com/KnutAM/FerriteAssembly.jl) 
is make it convenient to iterate over a grid and perform some task for each item.
The primary use is to iterate over all cells and assemble into global matrices and 
vectors. However, it makes many other tasks easy as well, for example iterating over faces for adding boundary conditions, or integrating a function over a given domain. 

The design is built around two main types of objects
1. **Workers**: What to do for a given item
2. **Domains**: What items to work on

The package user is responsible for writing the code for the actual work (e.g., writing the element routine). A given combination of a worker and a type of domain requires the user to overload a specific function. For example, an assembler (worker) and a cell domain requires overloading [`element_routine!`](@ref) or [`element_residual!`](@ref).

## Documentation structure

The documentation has two main parts

* Learning by doing: *Runnable examples*
  * Tutorials: Complete examples
  * Short guides for specific tasks
* API documentation: *Docstring and short explanations*
  * Domains: *What items to work on*
  * Workers: *What to do for an item*
  * Convenience: *Premade solutions to make life easy*
  * Customizations: *When the standard isn't good enough*
  * Internals: *Developer documentation*

## Heat equation example
```@eval
# Include the example here, but modify the Literate output to suit being embedded
using Literate, Markdown
filename = "firstexample_literate"
Literate.markdown(filename*".jl"; execute=true)
contents = read(filename*".md", String)
Literate.script(filename*".jl"; name="firstexample")
rm(filename*".jl")
rm(filename*".md")
header_end = last(findnext("```", contents, 4))+1
Markdown.parse(replace(contents[header_end:end], 
    "*This page was generated using [Literate.jl]"=>"*The examples were generated using [Literate.jl]")
    )
```
