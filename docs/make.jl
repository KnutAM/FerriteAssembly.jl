using FerriteAssembly
import FerriteAssembly: ExampleElements
import MaterialModelsBase as MMB
using Ferrite
using Test
using Documenter

const is_ci = get(ENV, "CI", "false") == "true"

include("generate.jl")
tutorials = ["heat_equation.jl", "viscoelasticity.jl", "incompressible_elasticity.jl", "mixed_materials.jl", "iga.jl"]
generated_tutorials = build_examples(tutorials; type="tutorials")
howto = [
    "threaded_assembly.jl", "automatic_differentiation.jl", "local_constraints.jl",
    "robin_bc.jl", "volume_integral.jl", "surface_integral.jl"]
generated_howto = build_examples(howto; type="howto")

DocMeta.setdocmeta!(FerriteAssembly, :DocTestSetup, :(using FerriteAssembly); recursive=true)

# Run example from `index.md` to force CI failure if it doesn't work
# include("src/firstexample_literate.jl")


makedocs(;
    modules=[FerriteAssembly],
    authors="Knut Andreas Meyer and contributors",
    repo=Remotes.GitHub("KnutAM", "FerriteAssembly.jl"),
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
        "Reference" => [
            "Design" => "design.md",
            "DomainBuffers" => [
                "DomainBuffers" => "DomainBuffers/Setup.md",
                "State variables" => "DomainBuffers/StateVariables.md",
                "ItemBuffer" => "DomainBuffers/ItemBuffer.md",
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
