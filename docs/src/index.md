```@meta
CurrentModule = FerriteAssembly
```

# FerriteAssembly
The goal of [FerriteAssembly](https://github.com/KnutAM/FerriteAssembly.jl) 
is to provide a simple structure for assembling in 
[Ferrite.jl](https://github.com/Ferrite-FEM/Ferrite.jl/).

## Key features
* Easy switching between *sequential* and *threaded* assembly
* *Multiple domains* with different fields, interpolations, and/or element routines.
* Efficient *automatic differentiation* if analytical tangent is not implemented. 
* Support for handling of (old and new) *state variables*
* Easy [*integration*](@ref Integration) of a function over the domain

## Typical workflow
1. Define your custom type and associated element routine (See [Example elements](@ref))
2. Setup Ferrite's `DofHandler` and `CellValues` as usual, and call [`setup_assembly`](@ref)
3. For each assembly, call [`doassemble!`](@ref)

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
