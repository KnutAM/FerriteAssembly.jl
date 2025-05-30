using FerriteAssembly
import FerriteAssembly: ExampleElements
import MaterialModelsBase as MMB
using Ferrite
using Test
using Documenter

const is_ci = get(ENV, "CI", "false") == "true"

include("download_assets.jl")

include("generate.jl")
tutorials = [
    "heat_equation.jl", 
    "viscoelasticity.jl", 
    "incompressible_elasticity.jl", 
    "mixed_materials.jl", 
    "iga.jl",
    "phasefield_fracture.jl"
    ]
generated_tutorials = build_examples(tutorials; type="tutorials")
howto = [
    "threaded_assembly.jl", "automatic_differentiation.jl", "local_constraints.jl",
    "robin_bc.jl", "volume_integral.jl", 
    "surface_integral.jl"]
generated_howto = build_examples(howto; type="howto")

DocMeta.setdocmeta!(FerriteAssembly, :DocTestSetup, :(using FerriteAssembly); recursive=true)

makedocs(;
    authors="Knut Andreas Meyer and contributors",
    sitename="FerriteAssembly.jl",
    format=Documenter.HTML(;
        canonical="https://KnutAM.github.io/FerriteAssembly.jl",
        assets=String[],
        collapselevel = 1,
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
                "Quadrature point eval" => "Workers/QuadPointEvaluator.md",
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
)

# Remove output files from build directory
clean_output_files(joinpath(@__DIR__, "build", "tutorials"))
clean_output_files(joinpath(@__DIR__, "build", "howto"))

deploydocs(;
    repo="github.com/KnutAM/FerriteAssembly.jl.git",
    devbranch="main",
    push_preview=true,
    versions = ["dev" => "dev", "stable" => "v^", "v#.#"] # dev first makes this default redirect
)
