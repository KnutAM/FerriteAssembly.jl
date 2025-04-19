

# Loading required packages 
using Ferrite, FerriteAssembly, Tensors, LinearAlgebra
using neper_jll, FerriteGmsh, DelimitedFiles
using MechanicalMaterialModels

#= 
We start by obtaining the Ferrite logo from neper, meshed with 
`Triangle` elements along with the grain orientation described 
as rotation vectors whose magnitude describes the rotation angle.
=#
grid, rotation = mktempdir() do tmp
    cd(tmp) do
        ## Tessalation to create the Ferite logo
        success(`$(neper()) -T -dim 2 -n 6 -id 4 -reg 1 -oridescriptor axis-angle:passive -o neperfile -format tess,ori`)
        ## Mesh the Ferrite logo
        meshed = success(`$(neper()) -M neperfile.tess -dim 2 -order 1 -rcl 2`)
        for i in 1:10 # Sometimes meshing fails, avoid CI failures due to this...   #hide
            meshed && break                                                         #hide
            @info "Retry nr $i/10"                                                  #hide
            meshed = success(`$(neper()) -M neperfile.tess -dim 2 -order 1 -rcl 2`) #hide
        end                                                                         #hide
        ## Read the mesh using FerriteGmsh
        grid = togrid("neperfile.msh")
        ## Clear the facet sets and add the boundaries
        empty!(grid.facetsets)
        addfacetset!(grid, "left", x -> abs(x[1]) < 1e-8)
        addfacetset!(grid, "bottom", x -> abs(x[2]) < 1e-8)
        addfacetset!(grid, "right", x -> x[1] ≈ 1)
        addfacetset!(grid, "top", x -> x[2] ≈ 1)
        ## Read the orientation file
        orientation = readdlm("neperfile.ori")
        rotation_directions = [Vec{3}(row[1:3]) for row in eachrow(orientation)]
        rotation_angles = deg2rad.(orientation[:, 4])
        rotation_vectors = [normalize(v) * α for (v, α) in zip(rotation_directions, rotation_angles)]
        ## "Return grid and orientation
        (grid, rotation_vectors)
    end
end

# Create the standard Ferrite setup 
ip = Lagrange{RefTriangle, 2}()^2
qr = QuadratureRule{RefTriangle}(2)
cv = CellValues(qr, ip)

dh = DofHandler(grid)
add!(dh, :u, ip)
close!(dh)

ch = ConstraintHandler(dh)
add!(ch, Dirichlet(:u, getfacetset(grid, "left"), Returns(0.0), 1)) # Fix x-displacement
add!(ch, Dirichlet(:u, getfacetset(grid, "bottom"), Returns(0.0), 2)) # Fix y-displacement
add!(ch, Dirichlet(:u, getfacetset(grid, "right"), (x, t) -> t, 1)) # Ramp uₓ proportional with time
close!(ch)

