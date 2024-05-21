using Ferrite, FerriteAssembly, BenchmarkTools
import FerriteAssembly.ExampleElements: StationaryFourier

struct StationaryFourierAD
    k::Float64 # Thermal conductivity
end
material = StationaryFourierAD(1.0)

function FerriteAssembly.element_residual!(re, state, ae,
        material::StationaryFourierAD, cellvalues, cellbuffer
        )
    n_basefuncs = getnbasefunctions(cellvalues)
    # Loop over quadrature points
    for q_point in 1:getnquadpoints(cellvalues)
        # Get the quadrature weight
        dΩ = getdetJdV(cellvalues, q_point)
        ∇u = function_gradient(cellvalues, q_point, ae)
        # Loop over test shape functions
        for i in 1:n_basefuncs
            ∇δN = shape_gradient(cellvalues, q_point, i)
            # re = fint - fext
            re[i] += material.k*(∇δN ⋅ ∇u) * dΩ
        end
    end
end;

grid = generate_grid(Quadrilateral, (100, 100))
ip = Lagrange{RefQuadrilateral,1}()
dh = DofHandler(grid); add!(dh, :u, ip); close!(dh)
cellvalues = CellValues(QuadratureRule{RefQuadrilateral}(2), ip);
K = create_sparsity_pattern(dh)
r = zeros(ndofs(dh));

domain = DomainSpec(dh, material, cellvalues)
buffer = setup_domainbuffer(domain);

a = zeros(ndofs(dh))
assembler = start_assemble(K, r)
work!(assembler, buffer; a=a);

domain_builtin = DomainSpec(dh, StationaryFourier(1.0), cellvalues)
buffer_builtin = setup_domainbuffer(domain_builtin);

if get(ENV, "CI", "false") == "true"            #hide
@btime work!($assembler, $buffer; a=$a)
@btime work!($assembler, $buffer_builtin; a=$a)
end                                             #hide

buffer_ad = setup_domainbuffer(domain; autodiffbuffer=true)
if get(ENV, "CI", "false") == "true"            #hide
@btime work!($assembler, $buffer_ad; a=$a)
end                                             #hide

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

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
