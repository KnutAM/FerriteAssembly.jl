using Ferrite

function get_random_cellnums(num_cells, split_fraction)
    num_split = round(Int, num_cells*split_fraction)
    @assert num_cells > num_split
    split_cells = sizehint!(Set{Int}(), num_split)
    while length(split_cells) < num_split
        push!(split_cells, rand(1:num_cells))
    end
    return sort!(collect(split_cells))
end

function generate_mixed_grid(nels::NTuple{2,Int}; split_fraction=0.5)
    g0 = generate_grid(Quadrilateral, nels)
    split_cells = get_random_cellnums(getncells(g0), split_fraction)
    new_cells = Union{Quadrilateral,Triangle}[cell for cell in getcells(g0)]
    old_cellset = Set(1:getncells(g0))
    new_cellset = sizehint!(Set{Int}(), 2*length(split_cells))
    for cellnr in split_cells
        old_cell = new_cells[cellnr]
        cell1 = Triangle(old_cell.nodes[1:3])
        cell2 = Triangle(old_cell.nodes[[1,3,4]])
        new_cells[cellnr] = cell1
        push!(new_cells, cell2)
        pop!(old_cellset, cellnr)
        push!(new_cellset, cellnr)
        push!(new_cellset, length(new_cells))
    end
    return Grid(new_cells, getnodes(g0), cellsets=Dict("old"=>old_cellset, "new"=>new_cellset))
end

function add_inclusion!(grid; size=0.5)
    addcellset!(grid, "inclusion", x -> all(v-> abs(v)<(size), x))
    addcellset!(grid, "matrix", setdiff(1:getncells(grid), getcellset(grid, "inclusion")))
end

function create_partitioned_grid()
    grid = generate_mixed_grid((20,20));
    add_inclusion!(grid)
end

function save_as_vtk()
    # TODO: Use FerriteViz instead and save in high quality pdf with CairoMakie
    grid = create_partitioned_grid()
    cell_type_data = zeros(Int, getncells(grid))
    cell_physical_domains = zeros(Int, getncells(grid))
    cell_work_domains = zeros(Int, getncells(grid))
    for (i, location) in enumerate(("matrix", "inclusion"))
        for k in getcellset(grid, location)
            cell_physical_domains[k] = i
        end
        for (j, cell_type) in enumerate(("old", "new"))
            for k in intersect(getcellset(grid, location), getcellset(grid, cell_type))
                cell_work_domains[k] = 2*(i-1) + j 
            end
            if i==1
                for k in getcellset(grid, cell_type)
                    cell_type_data[k] = j
                end
            end
        end
    end

    vtk_grid(joinpath(@__DIR__, "..", "..", "griddomains"), grid) do vtk
        vtk_cell_data(vtk, cell_type_data, "cell_type_data")
        vtk_cell_data(vtk, cell_physical_domains, "cell_physical_domains")
        vtk_cell_data(vtk, cell_work_domains, "cell_work_domains")
    end
end