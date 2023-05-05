module ExampleElements

using Ferrite    
import ..FerriteAssembly
    
include("WeakForm.jl")
include("HeatEquation.jl")
include("LinearElasticity.jl")
include("PorousMedia.jl")

end