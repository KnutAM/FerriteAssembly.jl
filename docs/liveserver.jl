#!/usr/bin/env julia

# Build the documentation and serve it locally, rebuilding whenever a source file
# changes. Pass `draft` to skip executing the tutorials and how-tos.

if !(isempty(ARGS) || ARGS == ["draft"])
    throw(ArgumentError("usage: julia docs/liveserver.jl [draft]"))
end

# Root of the repository
const repo_root = dirname(@__DIR__)

# Make sure docs environment is active and instantiated
import Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

# Communicate with docs/make.jl that it should serve without deploying
push!(ARGS, "liveserver")

# Run LiveServer.servedocs(...)
import LiveServer
LiveServer.servedocs(;
    host = "0.0.0.0",
    # Documentation root where make.jl and src/ are located
    foldername = joinpath(repo_root, "docs"),
    # Extra source folder to watch for changes
    include_dirs = [
        # Watch the src folder so docstrings can be Revise'd
        joinpath(repo_root, "src"),
    ],
    skip_dirs = [
        # Skip the folders where Literate.jl output is written. This is needed
        # to avoid infinite loops where running make.jl updates watched files,
        # which then triggers a new run of make.jl etc.
        joinpath(repo_root, "docs/src/tutorials"),
        joinpath(repo_root, "docs/src/howto"),
        # Skip the folder with downloaded assets, see docs/download_assets.jl
        joinpath(repo_root, "docs/src/assets"),
    ],
    include_files = [
        joinpath(repo_root, "docs/generate.jl"),
        joinpath(repo_root, "docs/download_assets.jl"),
    ],
)
