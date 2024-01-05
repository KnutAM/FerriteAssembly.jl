# generate examples
import Literate

function clean_output_files(dir)
    num_deleted = 0
    cd(dir) do 
        for file in readdir()
            if any(ext -> endswith(file, ext), (".vtu", ".pvd"))
                rm(file)
                num_deleted += 1
            end
        end
    end
    @info "Deleted $num_deleted output files from $dir"
end

# type="tutorials" or "howto" 
function build_examples(jl_files; type)
    EXAMPLEDIR = joinpath(@__DIR__, "src", "literate_$type")
    GENERATEDDIR = joinpath(@__DIR__, "src", type)
    rm(GENERATEDDIR; force=true, recursive=true)
    mkpath(GENERATEDDIR)

    # Copy supplementary files first
    supplementary_fileextensions = [".inp", ".svg", ".png", ".jpg", ".gif"]
    for example in readdir(EXAMPLEDIR)
        supplementary_ext = any(endswith.(example, supplementary_fileextensions))
        supplementary_jl = endswith(example, ".jl") && example ∉ jl_files
        if supplementary_ext || supplementary_jl
            cp(joinpath(EXAMPLEDIR, example), joinpath(GENERATEDDIR, example); force=true)
        end
    end

    for example in jl_files
        input = abspath(joinpath(EXAMPLEDIR, example))
        isfile(input) || throw(SystemError("$input not found"))
        script = Literate.script(input, GENERATEDDIR)
        code = strip(read(script, String))

        # remove "hidden" lines which are not shown in the markdown
        line_ending_symbol = occursin(code, "\r\n") ? "\r\n" : "\n"
        code_clean = join(filter(x->!endswith(x,"#hide"),split(code, r"\n|\r\n")), line_ending_symbol)

        mdpost(str) = replace(str, "@__CODE__" => code_clean)
        Literate.markdown(input, GENERATEDDIR, postprocess = mdpost)
        Literate.notebook(input, GENERATEDDIR, execute = is_ci) # Don't execute locally
    end
    clean_output_files(GENERATEDDIR)
    
    return map(f->joinpath(type, replace(f, ".jl"=>".md")), jl_files)
end
