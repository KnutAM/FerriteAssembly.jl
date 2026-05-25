# # Automatic differentiation
# **Executive summary:** Define `element_residual!` and call `setup_domainbuffer(domain; autodiffbuffer=true)`.
# 
# Here, we show how to define only the `element_residual!` function, which 
# does not calculate the element stiffness matrix, `Ke`, and then let `FerriteAssembly`'s
# Automatic Differentiation (AD) functionality (that builds on `ForwardDiff.jl`)
# calculate `Ke`. To compare with, we'll use the `StationaryFourier` example element that 
# does define an `element_routine!` that calculates `Ke` explicitly. 
using Ferrite, FerriteAssembly, BenchmarkTools
import FerriteAssembly.ExampleElements: StationaryFourier

# First, we'll define a new material, for which we only implement the `element_residual!`
struct StationaryFourierAD
    k::Float64 # Thermal conductivity
end
material = StationaryFourierAD(1.0)

function FerriteAssembly.element_residual!(re, state, ae, 
        material::StationaryFourierAD, cellvalues, cellbuffer
        )
    n_basefuncs = getnbasefunctions(cellvalues)
    ## Loop over quadrature points
    for q_point in 1:getnquadpoints(cellvalues)
        ## Get the quadrature weight
        dΩ = getdetJdV(cellvalues, q_point)
        ∇u = function_gradient(cellvalues, q_point, ae)
        ## Loop over test shape functions
        for i in 1:n_basefuncs
            ∇δN = shape_gradient(cellvalues, q_point, i)
            ## re = fint - fext
            re[i] += material.k*(∇δN ⋅ ∇u) * dΩ
        end
    end
end;

# Then, we create a standard Ferrite setup:
grid = generate_grid(Quadrilateral, (100, 100))
ip = Lagrange{RefQuadrilateral,1}()
dh = DofHandler(grid); add!(dh, :u, ip); close!(dh)
cellvalues = CellValues(QuadratureRule{RefQuadrilateral}(2), ip);
K = allocate_matrix(dh)
r = zeros(ndofs(dh));

# The assembly can then be done as usual
domain = DomainSpec(dh, material, cellvalues)
buffer = setup_domainbuffer(domain);

# But for doing automatic differentiation, we need 
# actual values for the degrees of freedom, so they 
# must be passed to the work function,
a = zeros(ndofs(dh))
assembler = start_assemble(K, r)
work!(assembler, buffer; a=a);

# ## Performance
# Using the simplest setup for AD is quite slow, due to lack of pre-allocations. 
# If we would use the builtin material instead, which doesn't use AD, it is much faster:
domain_builtin = DomainSpec(dh, StationaryFourier(1.0), cellvalues)
buffer_builtin = setup_domainbuffer(domain_builtin);
# Explicitly defining the element stiffness was a lot faster and has less allocations:
if get(ENV, "CI", "false") == "true"            #hide
@btime work!($assembler, $buffer; a=$a)
@btime work!($assembler, $buffer_builtin; a=$a)
end                                             #hide

# However, FerriteAssembly comes with a special `cellbuffer::AutoDiffCellBuffer` 
# for speeding up automatic differentiation. 
# Simply say that you want to use this during setup to improve the performance when 
# using AD
buffer_ad = setup_domainbuffer(domain; autodiffbuffer=true)
if get(ENV, "CI", "false") == "true"            #hide
@btime work!($assembler, $buffer_ad; a=$a)
end                                             #hide

# Test that all methods give the same stiffness #src
assembler = start_assemble(K, r)                #hide
work!(assembler, buffer; a=a);                  #hide
K1 = copy(K);                                   #hide
assembler = start_assemble(K, r)                #hide
work!(assembler, buffer_builtin; a=a);          #hide
K2 = copy(K);                                   #hide
assembler = start_assemble(K, r)                #hide
work!(assembler, buffer_ad; a=a);               #hide
K3 = copy(K);                                   #hide
                                                #hide
using Test                                      #hide
@test K1.nzval !== K2.nzval                     #hide
@test K1 ≈ K2 ≈ K3                              #hide
nothing;                                        #hide