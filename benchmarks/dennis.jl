using Ferrite, SparseArrays

using ThreadPinning
pinthreads(:cores)

function create_colored_cantilever_grid(celltype, n)
    grid = generate_grid(celltype, (10*n, n, n), Vec{3}((0.0, 0.0, 0.0)), Vec{3}((10.0, 1.0, 1.0)))
    colors = create_coloring(grid)
    return grid, colors
end;

# #### DofHandler
function create_dofhandler(grid::Grid{dim}) where {dim}
    dh = DofHandler(grid)
    push!(dh, :u, dim) # Add a displacement field
    close!(dh)
end;

# ### Stiffness tensor for linear elasticity
function create_stiffness(::Val{dim}) where {dim}
    E = 200e9
    ν = 0.3
    λ = E*ν / ((1+ν) * (1 - 2ν))
    μ = E / (2(1+ν))
    δ(i,j) = i == j ? 1.0 : 0.0
    g(i,j,k,l) = λ*δ(i,j)*δ(k,l) + μ*(δ(i,k)*δ(j,l) + δ(i,l)*δ(j,k))
    C = SymmetricTensor{4, dim}(g);
    return C
end;

# ## Threaded data structures
#
# ScratchValues is a thread-local collection of data that each thread needs to own,
# since we need to be able to mutate the data in the threads independently
struct ScratchValues{T, CV <: CellValues, FV <: FaceValues, TT <: AbstractTensor, dim, Ti}
    Ke::Matrix{T}
    fe::Vector{T}
    cellvalues::CV
    facevalues::FV
    global_dofs::Vector{Int}
    ɛ::Vector{TT}
    coordinates::Vector{Vec{dim, T}}
    assembler::Ferrite.AssemblerSparsityPattern{T, Ti}
end;

# Each thread need its own CellValues and FaceValues (although, for this example we don't use
# the FaceValues)
function create_values(refshape, dim, order::Int)
    ## Interpolations and values
    interpolation_space = Lagrange{refshape, 1}()
    quadrature_rule = QuadratureRule{refshape}(order)
    face_quadrature_rule = FaceQuadratureRule{refshape}(order)
    cellvalues = [CellValues(quadrature_rule, interpolation_space^dim) for i in 1:Threads.nthreads()];
    facevalues = [FaceValues(face_quadrature_rule, interpolation_space^dim) for i in 1:Threads.nthreads()];
    return cellvalues, facevalues
end;

# Create a `ScratchValues` for each thread with the thread local data
function create_scratchvalues(K, f, dh::DofHandler{dim}) where {dim}
    nthreads = Threads.nthreads()
    assemblers = [start_assemble(K, f) for i in 1:nthreads]
    cellvalues, facevalues = create_values(RefHexahedron, dim, 2)

    n_basefuncs = getnbasefunctions(cellvalues[1])
    global_dofs = [zeros(Int, ndofs_per_cell(dh)) for i in 1:nthreads]

    fes = [zeros(n_basefuncs) for i in 1:nthreads] # Local force vector
    Kes = [zeros(n_basefuncs, n_basefuncs) for i in 1:nthreads]

    ɛs = [[zero(SymmetricTensor{2, dim}) for i in 1:n_basefuncs] for i in 1:nthreads]

    coordinates = [[zero(Vec{dim}) for i in 1:length(dh.grid.cells[1].nodes)] for i in 1:nthreads]

    return [ScratchValues(Kes[i], fes[i], cellvalues[i], facevalues[i], global_dofs[i],
                         ɛs[i], coordinates[i], assemblers[i]) for i in 1:nthreads]
end;

# ## Threaded assemble

# The assembly function loops over each color and does a threaded assembly for that color
function doassemble(K::SparseMatrixCSC, colors, grid::Grid, dh::DofHandler, C::SymmetricTensor{4, dim}, f::Vector{Float64}, scratches::Vector{SV}, b::Vec{dim}) where {dim, SV}
    ## Each color is safe to assemble threaded
    for color in colors
        ## We try to equipartition the array to increase load per task.
        chunk_size = max(1, 1 + length(color) ÷ Threads.nthreads())
        color_partition = [color[i:min(i + chunk_size - 1, end)] for i in 1:chunk_size:length(color)]
        ## Now we should have a 1:1 correspondence between tasks and elements to assemble.
        Threads.@threads :static for i in 1:length(color_partition)
            for cellid ∈ color_partition[i]
                assemble_cell!(scratches[i], cellid, K, grid, dh, C, b)
            end
        end
    end

    return K, f
end

# The cell assembly function is written the same way as if it was a single threaded example.
# The only difference is that we unpack the variables from our `scratch`.
function assemble_cell!(scratch::ScratchValues, cell::Int, K::SparseMatrixCSC,
                        grid::Grid, dh::DofHandler, C::SymmetricTensor{4, dim}, b::Vec{dim}) where {dim}

    ## Unpack our stuff from the scratch
    Ke, fe, cellvalues, facevalues, global_dofs, ɛ, coordinates, assembler =
         scratch.Ke, scratch.fe, scratch.cellvalues, scratch.facevalues,
         scratch.global_dofs, scratch.ɛ, scratch.coordinates, scratch.assembler

    fill!(Ke, 0)
    fill!(fe, 0)

    n_basefuncs = getnbasefunctions(cellvalues)

    ## Fill up the coordinates
    nodeids = grid.cells[cell].nodes
    for j in 1:length(coordinates)
        coordinates[j] = grid.nodes[nodeids[j]].x
    end

    reinit!(cellvalues, coordinates)

    for q_point in 1:getnquadpoints(cellvalues)
        for i in 1:n_basefuncs
            ɛ[i] = symmetric(shape_gradient(cellvalues, q_point, i))
        end
        dΩ = getdetJdV(cellvalues, q_point)
        for i in 1:n_basefuncs
            δu = shape_value(cellvalues, q_point, i)
            fe[i] += (δu ⋅ b) * dΩ
            ɛC = ɛ[i] ⊡ C
            for j in 1:n_basefuncs
                Ke[i, j] += (ɛC ⊡ ɛ[j]) * dΩ
            end
        end
    end

    celldofs!(global_dofs, dh, cell)
    assemble!(assembler, global_dofs, fe, Ke)
end;

function setup_assemble()
    refshape = RefHexahedron
    dim = 3
    n = 30
    grid, colors = create_colored_cantilever_grid(Hexahedron, n);
    dh = create_dofhandler(grid);

    K = create_sparsity_pattern(dh);
    C = create_stiffness(Val{3}());
    f = zeros(ndofs(dh))
    scratches = create_scratchvalues(K, f, dh)
    b = Vec{3}((0.0, 0.0, 0.0))
    ## compilation
    return K, colors, grid, dh, C, f, scratches, b
end

function timeit(args...)
    return @elapsed doassemble(args...)
end

args = setup_assemble()
doassemble(args...) # compilation

open(joinpath(@__DIR__, "dennis.txt"); append=true) do io
    write(io, string(Threads.nthreads()))
    for _ in 1:3
        t = timeit(args...)
        write(io, ", ", string(t))
    end
    write(io, "\n")
end

