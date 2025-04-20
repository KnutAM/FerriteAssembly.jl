using Downloads: download

function download_assets(storepath::String, assets::Vector{String})
    mkpath(storepath)
    for asset in assets
        afile = joinpath(storepath, asset)
        isfile(afile) || download(FerriteAssembly.asset_url(asset), afile)
    end
end

download_assets(
    joinpath(@__DIR__, "src", "assets"),
    ["domains.svg"]
    )

download_assets(
    joinpath(@__DIR__, "src", "literate_tutorials"), 
    [
        "mixed_materials.png", "zener.svg", 
        "sent_results.png", "sent_animation.gif"
    ],
    )
