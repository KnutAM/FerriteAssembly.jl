```@meta
CurrentModule = FerriteAssembly
```

# FerriteAssembly
The goal of [FerriteAssembly](https://github.com/KnutAM/FerriteAssembly.jl) 
is make it convenient to iterate over a grid and perform some task for each item.
The primary use is to iterate over all cells and assemble into global matrices and 
vectors. However, it makes many other tasks easy as well, for example iterating over faces for adding boundary conditions, or integrating a function over a given domain. 

The design is built around two main object types
1. **Worker**: What to do for each item
2. **Domain**: What items to work on

More details are given in [Package Design](design.md), which is a recommended read for understanding the structure of the package. The current page highlights an example for doing standard assembly over cells. But before diving into that example, we give a short list with links to examples for some common tasks that `FerriteAssembly.jl` can help you with.

**TODO:** *Add links and examples*

* Easy Neumann boundary conditions
* Models with state variables
* Multiple domains
* Mixed grids
* Calculate integral quantities
* Postprocess data in each cell
* Robin boundary conditions

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
