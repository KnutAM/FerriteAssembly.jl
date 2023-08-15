using FerriteAssembly
import FerriteAssembly: ExampleElements
import MaterialModelsBase as MMB
using Test
using Documenter

const is_ci = get(ENV, "CI", "false") == "true"

include("generate.jl")
tutorials = ["plasticity.jl", "mixed_materials.jl"]
generated_tutorials = build_examples(tutorials; type="tutorials")
howto = ["robin_bc.jl", "volume_integral.jl", "surface_integral.jl"]
generated_howto = build_examples(howto; type="howto")

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
        "Learning by doing" => [
            "Tutorials" => generated_tutorials,
             "How-to" => generated_howto
        ],
        "API documentation" => [
            "Design" => "design.md",
            "DomainBuffers" => [
                "DomainBuffers" => "DomainBuffers/Setup.md",
                "State variables" => "DomainBuffers/StateVariables.md",
                "CellBuffer" => "DomainBuffers/CellBuffer.md",
            ],
            "Workers" => [
                "Workers" => "Workers/Workers.md",
                "Assemblers" => "Workers/Assemblers.md",
                "Integrators" => "Workers/Integrators.md",
            ],
            "Convenience" => [
                "External loads" => "Convenience/LoadHandler.md",
                "Mechanical materials" => "Convenience/MaterialModelsBase.md",
                "Example elements" => "Convenience/ExampleElements.md",
            ],
            "Customizations" => "Customization.md",
            "Internals" => "internals.md",
        ]
    ],
    #strict=true,
)

deploydocs(;
    repo="github.com/KnutAM/FerriteAssembly.jl",
    devbranch="main",
    push_preview=true,
)
