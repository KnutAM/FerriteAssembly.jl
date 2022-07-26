using FerriteAssembly
using Documenter

const is_ci = get(ENV, "CI", "false") == "true"

include("generate.jl")
examples = ["plasticity.jl", "mixed_materials.jl"]
GENERATEDEXAMPLES = [joinpath("examples", replace(f, ".jl"=>".md")) for f in examples]

build_examples(examples)

DocMeta.setdocmeta!(FerriteAssembly, :DocTestSetup, :(using FerriteAssembly); recursive=true)


# Run example from `index.md` to force CI failure if it doesn't work
include("src/firstexample_literate.jl")


makedocs(;
    modules=[FerriteAssembly],
    authors="Knut Andreas Meyer and contributors",
    repo="https://github.com/KnutAM/FerriteAssembly.jl/blob/{commit}{path}#{line}",
    sitename="FerriteAssembly.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://KnutAM.github.io/FerriteAssembly.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Examples" => GENERATEDEXAMPLES,
        "API" => "api.md",
        "Datastructures" => "datastructures.md",
        "Internals" => "internals.md",
    ],
    strict=true,
)

deploydocs(;
    repo="github.com/KnutAM/FerriteAssembly.jl",
    devbranch="main",
    push_preview=true,
)
