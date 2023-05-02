module ExampleElements

using Ferrite    
import ..FerriteAssembly as FA
    
include("HeatEquation.jl")
include("LinearElasticity.jl")
include("PorousMedia.jl")

end