"""
    DomainSpec(sdh::SubDofHandler, material, fe_values; set=_getcellset(sdh), colors_or_chunks=nothing, user_data=nothing)
    DomainSpec(dh::DofHandler, material, fe_values; set=1:getncells(dh), colors=nothing, chunks=nothing, user_data=nothing)

Create a `DomainSpec` that can be used to set up a domain buffer. 

* `sdh`/`dh`: Give the `DofHandler` for the domain in question, or a `SubDofHandler` in case there are more than one in `DofHandler` 
  (See `Ferrite.jl's documentation`)
* `material`: Used for dispatch on the utilized `worker`'s function. 
* `fe_values`: `CellValues` or `FacetValues` depending on the type of domain. 
* `set`: The items in the domain, the element type determines the type of domain 
  - Cell domain: `Int`
  - Facet domain: `FacetIndex`
* `colors::Vector{Vector{I}}`: used to avoid race conditions when multithreading. 
  For cell domains, `I=Int`, and for facet domains, `I` can be either `Int` 
  (denoting cell numbers) or `FacetIndex` for actual facets. If `I=Int`, it will be 
  converted to `FacetIndex` internally. If `colors=nothing` and `chunks=nothing`, 
  `Ferrite.jl`'s default coloring algorithm is used.
* `chunks::Vector{Vector{Vector{I}}}`. During multithreading, each task works with 
  items in one `chunk::Vector{I}` at a time. Items in `chunks[k][i]` and `chunks[k][j]`
  should be independent (i.e. not share dofs). If given, this input takes precedence over 
  `colors`. For `chunks`, `I` must be `Int` for cell domains and `FacetIndex` for facet domains. 
* `user_data`: Can be whatever the user wants to and is passed along by reference everywhere. 
  It is accessible from the `ItemBuffer` (e.g. `CellBuffer`) given to the `worker`'s function
  via the `get_user_data` function. However, since it is passed by reference, modifying values 
  during work, care must be taken to ensure thread safety. 
  To avoid allocations, caches can be created separately with `allocate_cell_cache` and 
  `allocate_facet_cache`.
"""
struct DomainSpec{I}
    sdh::SubDofHandler
    material::Any
    fe_values::Any
    set::Union{AbstractVector{I},AbstractSet{I}}
    colors_or_chunks::Any
    user_data::Any
end
function DomainSpec(sdh::SubDofHandler, material, values; set=_getcellset(sdh), colors=nothing, chunks=nothing, user_data=nothing)
    intersected_set = intersect_cellset_sort(set,_getcellset(sdh))
    colors_or_chunks = chunks!==nothing ? chunks : colors
    return DomainSpec(sdh, material, values, intersected_set, colors_or_chunks, user_data)
end
function DomainSpec(dh::DofHandler, args...; kwargs...)
    length(dh.subdofhandlers) > 1 && error("A domain requires a single SubDofHandler")
    return DomainSpec(dh.subdofhandlers[1], args...; kwargs...)
end

"""
    setup_domainbuffers(domains::Dict{String,DomainSpec}, suppress_warnings = false; kwargs...)

Setup multiple domain buffers, one for each `DomainSpec` in `domains`.
Set `suppress_warnings = true` to suppress warnings checking for typical input errors when setting up multiple domains. 
See [`setup_domainbuffer`](@ref) for description of the keyword arguments.
"""
function setup_domainbuffers(domains::Dict{String,<:DomainSpec}, suppress_warnings = false; kwargs...)
    dbs = Dict(name => setup_domainbuffer(domain; kwargs...) for (name, domain) in domains)
    suppress_warnings || check_input(dbs)
    return dbs 
end

check_input(dbs::DomainBuffers) = check_input(dbs, eltype(getset(first(values(dbs)))))

function check_input(dbs::DomainBuffers, ::Type{I}) where {I <: Integer} # Domain of cells
    if !all(db -> eltype(getset(db)) === I, values(dbs))
        @warn "Not all sets have the same eltype, this is likely to cause errors - proceed with caution"
    end

    grid = get_grid(dbs)
    
    # Check that all cells are included, and
    num_cells_in_sets = sum(length ∘ getset, values(dbs); init = 0)
    num_cells_in_grid = getncells(grid)
    more_fewer = num_cells_in_sets > num_cells_in_grid ? "more" : "fewer"
    if num_cells_in_grid != num_cells_in_sets
        @warn "There are $more_fewer cells ($num_cells_in_sets) assigned to domains than cells in the grid ($num_cells_in_grid)"
    end
    included = zeros(Bool, getncells(grid))
    for db in values(dbs)
        for i in getset(db)
            included[i] = true
        end
    end
    num_missing = getncells(grid) - sum(included)
    all(included) || @warn "$num_missing cells are not included in the domainbuffers"
end

function check_input(dbs::DomainBuffers, ::Type) # Unspecified (typically facet set) without tests
    # TODO: Add check for non-disjoint sets (always applicable)
    return nothing 
end

"""
    setup_domainbuffer(domain::DomainSpec; a=nothing, threading=false, autodiffbuffer=false, num_tasks=Threads.nthreads())

Setup a domain buffer for a single grid domain, `domain`. 
* `a::Vector`: The global degree of freedom values are used to pass the 
  local element dof values to the [`create_cell_state`](@ref) function, 
  making it possible to create the initial state dependent on the initial 
  conditions for the field variables. 
* `threading`: Should a `ThreadedDomainBuffer` be created to work the grid multithreaded
  if supported by the used `worker`?
* `num_tasks`: The number of tasks to spawn during threaded assembly. Only applicable for `threading = true`.
* `autodiffbuffer`: Should a custom itembuffer be used to speed up the automatic 
  differentiation (if supported by the itembuffer)
"""
function setup_domainbuffer(domain::DomainSpec; threading=Val(false), kwargs...)
    return _setup_domainbuffer(threading, domain; kwargs...)
end

create_states(domain::DomainSpec{Int}, a) = create_states(domain.sdh, domain.material, domain.fe_values, a, domain.set, create_dofrange(domain.sdh))
create_states(::DomainSpec{FacetIndex}, ::Any) = Dict{Int,Nothing}()

function setup_itembuffer(adb, domain::DomainSpec{FacetIndex}, args...)
    dofrange = create_dofrange(domain.sdh)
    return setup_facetbuffer(adb, domain.sdh, domain.fe_values, domain.material, dofrange, domain.user_data)
end
function setup_itembuffer(adb, domain::DomainSpec{Int}, states)
    dofrange = create_dofrange(domain.sdh)
    return setup_cellbuffer(adb, domain.sdh, domain.fe_values, domain.material, first(values(states)), dofrange, domain.user_data)
end

function _setup_domainbuffer(threaded, domain; a=nothing, autodiffbuffer=Val(false), kwargs...)
    new_states = create_states(domain, a)
    old_states = create_states(domain, a)
    itembuffer = setup_itembuffer(autodiffbuffer, domain, new_states)
    return _setup_domainbuffer(threaded, domain.set, itembuffer, StateVariables(old_states, new_states), domain.sdh, domain.colors_or_chunks; kwargs...)
end

# Type-unstable switch
function _setup_domainbuffer(threaded::Bool, args...; kwargs...)
    return _setup_domainbuffer(Val(threaded), args...; kwargs...)
end
# Sequential
function _setup_domainbuffer(::Val{false}, set, itembuffer, states, sdh, args...; kwargs...)
    return DomainBuffer(set, itembuffer, states, sdh; kwargs...)
end
# Threaded
function _setup_domainbuffer(::Val{true}, set, itembuffer, states, sdh, args...; kwargs...)
    return ThreadedDomainBuffer(set, itembuffer, states, sdh, args...; kwargs...)
end