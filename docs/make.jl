using FerriteAssembly
using Documenter

DocMeta.setdocmeta!(FerriteAssembly, :DocTestSetup, :(using FerriteAssembly); recursive=true)

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
    ],
)

deploydocs(;
    repo="github.com/KnutAM/FerriteAssembly.jl",
    devbranch="main",
)
