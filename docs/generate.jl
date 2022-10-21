# generate examples
import Literate

function replace_include_function(filename::String)
    pat = "include(\"$filename\")"
    rep = read(joinpath(@__DIR__, "src", "literate", filename), String)
    return str -> replace(str, pat=>rep)
end


function build_examples(examples)
    EXAMPLEDIR = joinpath(@__DIR__, "src", "literate")
    GENERATEDDIR = joinpath(@__DIR__, "src", "examples")
    rm(GENERATEDDIR; force=true, recursive=true)
    mkpath(GENERATEDDIR)

    # Copy supplementary files first
    suplementary_fileextensions = [".inp", ".svg", ".png", ".jpg", ".gif"]
    for example in readdir(EXAMPLEDIR)
        supplementary_ext = any(endswith.(example, suplementary_fileextensions))
        supplementary_jl = endswith(example, ".jl") && example ∉ examples
        if supplementary_ext || supplementary_jl
            cp(joinpath(EXAMPLEDIR, example), joinpath(GENERATEDDIR, example); force=true)
        end
    end

    for example in examples
        input = abspath(joinpath(EXAMPLEDIR, example))
        isfile(input) || throw(SystemError("$input not found"))
        script = Literate.script(input, GENERATEDDIR)
        code = strip(read(script, String))

        # remove "hidden" lines which are not shown in the markdown
        line_ending_symbol = occursin(code, "\r\n") ? "\r\n" : "\n"
        code_clean = join(filter(x->!endswith(x,"#hide"),split(code, r"\n|\r\n")), line_ending_symbol)

        mdpost(str) = replace(str, "@__CODE__" => code_clean)
        prepost = example!="plasticity.jl" ? identity :
            replace_include_function("J2Plasticity.jl") ∘ replace_include_function("MaterialModelsBaseElement.jl")
        Literate.markdown(input, GENERATEDDIR, preprocess = prepost, postprocess = mdpost)
        Literate.notebook(input, GENERATEDDIR, execute = is_ci) # Don't execute locally
    end

    cd(GENERATEDDIR) do
        foreach(file -> any(endswith(file, ext) for ext in (".vtu", ".pvd", ".jld2")) && rm(file), readdir())
    end
end
