module ExampleElements

using Ferrite    
import ..FerriteAssembly
    
include("HeatEquation.jl")
include("LinearElasticity.jl")
include("PorousMedia.jl")

end