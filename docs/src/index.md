```@meta
CurrentModule = FerriteAssembly
```

# FerriteAssembly

The goal of [FerriteAssembly](https://github.com/KnutAM/FerriteAssembly.jl) 
is to provide a simple structure to perform assembly in 
[Ferrite.jl](https://github.com/Ferrite-FEM/Ferrite.jl/).

Sequential and threaded assembly of both the `DofHandler` and `MixedDofHandler` are supported.

The package works by exporting the `doassemble!` function, and require the 
user to define either `element_routine!` (calculate both `Ke` and `re`),
or just `element_residual!` (calculate only `re`). 
In the latter case, `Ke` is calculated by 
[ForwardDiff.jl](https://github.com/JuliaDiff/ForwardDiff.jl)

Dispatch is typically done on a user-defined `material` struct,
and possible as well on `cellvalues` (potentially a `NamedTuple/Tuple`).
State variables and current dof-values for the cell are directly available in the `element_routine!`
Old dof-values for the cell, user-defined `cache` and `cell_load` types, cell coordinates and more 
are available through the `CellBuffer` type given as additional input. 

## A minimal example
```@eval
# Include the example here, but modify the Literate output to suit being embedded
using Literate, Markdown
filename = "firstexample_literate"
Literate.markdown(filename*".jl")
contents = read(filename*".md", String)
Literate.script(filename*".jl"; name="firstexample")
rm(filename*".jl")
rm(filename*".md")
header_end = last(findnext("```", contents, 4))+1
Markdown.parse(replace(contents[header_end:end], 
    "*This page was generated using [Literate.jl]"=>"*The examples were generated using [Literate.jl]")
    )
```

## Detailed API description
### Overloaded element routine
One of the element methods should be overloaded for a given combination of `cellvalues`
and `material`. 
Note that the `cellvalues` are already `reinit!`:ed when passed to the element routines. 
```@docs
FerriteAssembly.element_routine!
FerriteAssembly.element_residual!
```

### `CellBuffer`
Variables that are used and modified for each cell of a certain type, 
but that don't belong to a specific cell, are collected in a `CellBuffer`.
```@docs
CellBuffer
```

### State variables
The initial state variables may vary depending on the position in the grid.
Furthermore, the datastructure depends on the type of dof handler, so
a convenience function exists that creates the correct variable
```@docs
create_states
```

### `doassemble!`
Once everything is set up, one can call the function which will actually 
do the assembly:
```@docs
doassemble!
```

### Threaded assembly
For parallel assembly, we need a vector of `CellBuffer`s and assemblers, 
one for each thread. 

For the `MixedDofHandler`, the outer loop is over the type of cells,
so we need a tuple that contains vectors of `CellBuffer`s. 
Construction of this via `deepcopy` is implemented as
```@docs
create_threaded_CellBuffers
```

A vector of assemblers that is convieniently created by calling 
```@docs
create_threaded_assemblers
```
