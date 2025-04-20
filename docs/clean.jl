
# Delete build directory
docpath(args...) = joinpath(@__DIR__, args...)
cleanpath(args...) = rm(joinpath(@__DIR__, args...); recursive = true, force = true)
cleanpath("build")
cleanpath.(("src",), ("tutorials", "howto", "assets"))
