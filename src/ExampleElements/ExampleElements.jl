module ExampleElements

using Ferrite
using MaterialModelsBase
import ..FerriteAssembly
    
include("WeakForm.jl")
include("HeatEquation.jl")
include("LinearElasticity.jl")
include("PorousMedia.jl")

include("MaterialModels/J2Plasticity.jl")

end