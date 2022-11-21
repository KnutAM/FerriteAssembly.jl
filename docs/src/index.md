```@meta
CurrentModule = FerriteAssembly
```

# FerriteAssembly

The goal of [FerriteAssembly](https://github.com/KnutAM/FerriteAssembly.jl) 
is to provide a simple structure for assembling in 
[Ferrite.jl](https://github.com/Ferrite-FEM/Ferrite.jl/).

Sequential and threaded assembly when using either the `DofHandler` or the `MixedDofHandler`,
including a possibility of mixed materials, is supported

The package works by exporting the [`doassemble!`](@ref) function, and requiring the 
user to define either [`element_routine!`](@ref) (calculate both `Ke` and `re`),
or just [`element_residual!`](@ref) (calculate only `re`). 
In the latter case, `Ke` is calculated by 
[ForwardDiff.jl](https://github.com/JuliaDiff/ForwardDiff.jl)

Dispatch `element_routine!`/`element_residual!` is typically based 
on a user-defined `material` struct, and possibly also on `cellvalues`.
For multiple fields, the latter can also be a `NamedTuple/Tuple` 
(or any other type that supports `Ferrite.reinit!`).
State variables (to be mutated) and current dof-values for the cell 
are directly available in the `element_routine!`

Old dof-values for the cell, user-defined `cache` and `cell_load` types, 
cell coordinates and more are available through the `CellBuffer` type that
is given as an additional input.

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
