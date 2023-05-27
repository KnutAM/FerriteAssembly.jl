include(joinpath(@__DIR__, "threaded_heatequation.jl"))
@time doassemble!(assembler, new_states, buffer)
@time doassemble!(assembler, new_states, buffer)
@time doassemble!(assembler, new_states, buffer)

