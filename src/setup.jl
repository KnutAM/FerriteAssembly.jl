"""
    DomainSpec(sdh::SubDofHandler, material, fe_values; set=_getcellset(sdh), colors=nothing, chunks=nothing, user_data=nothing, algebraic_variables=())
    DomainSpec(dh::DofHandler, material, fe_values; set=1:getncells(dh), colors=nothing, chunks=nothing, user_data=nothing, algebraic_variables=())

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
* `algebraic_variables::Tuple`: **Experimental**. The algebraic variables (added to the dofhandler
  with `add!(dh, name, AlgebraicVariable())`) used in the domain. Each element is either
  - `name::Symbol`: The degrees of freedom of the variable are added after the degrees of freedom
    of all fields for each item (cell or facet), see [Algebraic variables](@ref).
  - [`AttachAlgebraic(name => field)`](@ref FerriteAssembly.AttachAlgebraic): The degrees of freedom of
    the variable are added directly after the degrees of freedom of `field` for each item,
    see [Attached variables](@ref).
"""
struct DomainSpec{I}
    sdh::SubDofHandler
    material::Any
    fe_values::Any
    set::Union{AbstractVector{I},AbstractSet{I}}
    colors_or_chunks::Any
    user_data::Any
    algebraic_variables::Tuple
end
function DomainSpec(sdh::SubDofHandler, material, values; set=_getcellset(sdh), colors=nothing, chunks=nothing, user_data=nothing, algebraic_variables=())
    intersected_set = intersect_cellset_sort(set,_getcellset(sdh))
    colors_or_chunks = chunks!==nothing ? chunks : colors
    if !isa(algebraic_variables, Tuple{Vararg{Union{Symbol, AttachAlgebraic}}})
        throw(ArgumentError(string("`algebraic_variables` must be a tuple with the names (`Symbol`s) of the ",
            "algebraic variables or `AttachAlgebraic` declarations, got $(repr(algebraic_variables))")))
    end
    return DomainSpec(sdh, material, values, intersected_set, colors_or_chunks, user_data, algebraic_variables)
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

**Note:** All domains in `domains` must be of the same kind: either all cell domains or
all facet (or other non-cell) domains. Mixing cell and facet domains in a single
dictionary is not supported; use separate dictionaries (and separate calls to
`setup_domainbuffers`/`work!`) instead.
"""
function setup_domainbuffers(domains::Dict{String,<:DomainSpec}, suppress_warnings = false; kwargs...)
    isempty(domains) && return DomainBuffers{AbstractDomainBuffer}()
    if !allequal(_domain_entity, values(domains))
        throw(ArgumentError("Not allowed to mix different entities (e.g. cell and facet) in the same domainbuffer dictionary"))
    end
    dbs = DomainBuffers(name => setup_domainbuffer(domain; kwargs...) for (name, domain) in domains)
    suppress_warnings || check_input(dbs)
    return dbs
end

_domain_entity(::DomainSpec{I}) where {I} = _domain_entity(I)
_domain_entity(::DomainBuffer{I}) where {I} = _domain_entity(I)
_domain_entity(::ThreadedDomainBuffer{I}) where {I} = _domain_entity(I)

_domain_entity(::Type{<:Integer}) = :cell
_domain_entity(::Type{<:FacetIndex}) = :facet
# _domain_entity(::Type{<:InterfaceIndex}) # TODO, see #22

function check_input(dbs::DomainBuffers)
    isempty(dbs) && return nothing
    if all(db -> _domain_entity(db) === :cell, values(dbs))
        check_cell_domains(dbs)
    else
        check_other_domains(dbs)
    end
    return nothing
end

# Domains of cells: check eltype consistency and that all grid cells are covered exactly once.
function check_cell_domains(cell_dbs::DomainBuffers)
    grid = get_grid(cell_dbs)

    num_cells_in_sets = sum(length ∘ get_set, values(cell_dbs); init = 0)
    num_cells_in_grid = getncells(grid)
    more_fewer = num_cells_in_sets > num_cells_in_grid ? "more" : "fewer"
    if num_cells_in_grid != num_cells_in_sets
        @warn "There are $more_fewer cells ($num_cells_in_sets) assigned to domains than cells in the grid ($num_cells_in_grid)"
    end
    included = zeros(Bool, getncells(grid))
    for db in values(cell_dbs)
        for i in get_set(db)
            included[i] = true
        end
    end
    num_missing = getncells(grid) - sum(included)
    all(included) || @warn "$num_missing cells are not included in the domainbuffers"
end

# Domains of e.g. facets: no checks implemented yet.
function check_other_domains(other_dbs::DomainBuffers)
    @assert allequal(_domain_entity, values(other_dbs)) # Must be same, see `setup_domainbuffers`
    # (should never happen since error should be thrown earlier)
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
An empty `domain.set` (e.g. because the given `set` does not intersect the SubDofHandler's
cellset) is supported as a no-op: the resulting buffer is skipped when `work!`ing.
"""
function setup_domainbuffer(domain::DomainSpec; threading=Val(false), kwargs...)
    return _setup_domainbuffer(threading, domain; kwargs...)
end

DofBuffer(domain::DomainSpec) = DofBuffer(domain.sdh, setup_layout(domain.sdh, domain.algebraic_variables))

function create_states(domain::DomainSpec{Int}, dofbuffer::DofBuffer, a)
    return create_states(dofbuffer, domain.sdh, domain.material, domain.fe_values, a, domain.set)
end
create_states(::DomainSpec{FacetIndex}, ::Any, ::Any) = Dict{Int,Nothing}()

# Type-unstable switch
setup_itembuffer(autodiff::Bool, domain::DomainSpec, dofbuffer::DofBuffer, states) = setup_itembuffer(Val(autodiff), domain, dofbuffer, states)
function setup_itembuffer(::Val{false}, domain::DomainSpec{Int}, dofbuffer::DofBuffer, states)
    # An empty `states` (empty domain.set) has no real cell to sample; `nothing` is used as a
    # placeholder, matching facet buffers, since it is never dereferenced (work! is a no-op).
    sample_state = isempty(states) ? nothing : first(values(states))
    coords = sample_coordinates(domain.sdh)
    return CellBuffer(dofbuffer, coords, domain.fe_values, domain.material, sample_state, domain.user_data)
end
function setup_itembuffer(::Val{true}, domain::DomainSpec{Int}, dofbuffer::DofBuffer, states)
    return AutoDiffCellBuffer(setup_itembuffer(Val(false), domain, dofbuffer, states))
end
function setup_itembuffer(::Val{false}, domain::DomainSpec{FacetIndex}, dofbuffer::DofBuffer, states)
    coords = sample_coordinates(domain.sdh)
    return FacetBuffer(dofbuffer, coords, domain.fe_values, domain.material, domain.user_data)
end
function setup_itembuffer(::Val{true}, ::DomainSpec{FacetIndex}, ::DofBuffer, states)
    error("AutoDiffBuffer not implemented for FacetBuffer")
end

function _setup_domainbuffer(threaded, domain; a=nothing, autodiffbuffer=Val(false), kwargs...)
    dofbuffer = DofBuffer(domain)
    check_attached_values(dofbuffer.layout, domain.sdh, domain.fe_values)
    new_states = create_states(domain, dofbuffer, a)
    old_states = create_states(domain, dofbuffer, a)
    itembuffer = setup_itembuffer(autodiffbuffer, domain, dofbuffer, new_states)
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

# The values for fields with attached variables must have one shape function for each local dof
check_attached_values(::Nothing, ::SubDofHandler, ::Any) = nothing
function check_attached_values(layout::AlgebraicLayout, sdh::SubDofHandler, values)
    layout.fields === nothing && return nothing
    fieldnames = keys(layout.field_ranges)
    for (field, range) in pairs(layout.field_ranges)
        length(range) == length(dof_range(sdh, field)) && continue # No attached dofs
        field_values = if values isa NamedTuple
            haskey(values, field) || throw(ArgumentError(string("Variables are attached to :$field, but the values ",
                "(a NamedTuple with keys $(keys(values))) do not contain values for :$field")))
            values[field]
        elseif length(fieldnames) == 1
            values
        else
            throw(ArgumentError(string("Variables are attached to :$field, and the values for each field must then ",
                "be given as a NamedTuple with the field names as keys")))
        end
        n = getnbasefunctions(field_values)
        if n != length(range)
            throw(ArgumentError(string("The values for :$field have $n base functions, but the field has ",
                "$(length(range)) local dofs, including those of the attached variables")))
        end
    end
    return nothing
end
