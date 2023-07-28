"""
    Integrator(val::V; domains=nothing)

Integrate over the domain by modifying the value `val` for each cell in the function 
`integrate_cell`, which should be overloaded for the specific combination of `val::V` 
and the material for that cell. `V` must be a mutable type.

Specify `domains` to only integrate over a part of the grid.
It should contain the names of the `AssemblyDomains` to integrate over.
Single, `String`, or collections; `Set{String}`, `AbstractVector{String}`, or `NTuple{N,String}`
inputs are supported. If `domains=nothing`, all domains are integrated. 

**Note**: This integrator will currently only run sequentially even if a 
threaded assembly is setup with [`setup_assembly`](@ref)
"""
struct Integrator{T,D<:Union{Nothing,Set{String}}}
    val::T
    domains::D # Which domains to integrate over
end
Integrator(val; domains=nothing) = Integrator(val, domains)
Integrator(val, domains::Union{AbstractVector{String},NTuple{<:Any,String}}) = Integrator(val, Set(domains))
Integrator(val, domain::String) = Integrator(val, Set((domain,)))
# Properties:
can_thread(::Integrator) = false # Might support in future, but requires user-facing interface 
need_colors(::Integrator) = false
skip_this_domain(::Integrator{<:Any,Nothing}, ::String) = false
skip_this_domain(ig::Integrator{<:Any,<:Set}, name::String) = name ∉ ig.domains

"""
    integrate_cell!(val, cell_state, ae, material, cv, cellbuffer)

Function to be overloaded for `val` (potentially in combination with `material`),
and will be called when using [`Integrator`](@ref). Mutate `val` to add to the result.

The order of the inputs is chosen to follow the element routines
"""
function integrate_cell! end

function work_single_cell!(assembler::Integrator, cellbuffer)
    cv = get_values(cellbuffer)
    m = get_material(cellbuffer)
    cell_state = get_new_state(cellbuffer)
    ae = get_ae(cellbuffer)
    integrate_cell!(assembler.val, cell_state, ae, m, cv, cellbuffer)
end


"""
    SimpleIntegrator(fun::Function, val)

Calculate the integral
```math
\\int_\\Omega f(u, \\nabla u, s)\\, \\mathrm{d}\\Omega 
```
For single-field problems, the function signature of ``f`` is `fun(u, ∇u, qp_state)`.
`u` is the current function value, `∇u` the current function gradient, and `qp_state` the 
current state in the current quadrature point. This assumes that `cell_state::AbstractVector`, 
otherwise, `qp_state = cell_state`.  

For multi-field problem, we have `fun(u::NamedTuple, ∇u::NamedTuple, qp_state)`, 
where the keys in `u` and `∇u` are the fieldnames. The rest is same as for single-field 
problems. 

It is the user's responsibility that `fun(args...)::typeof(val)`. 
Additionally, the type of `val` must support
- `val + val`
- `val * x` (where `x::Real`)
- `zero(val)`
One exception to these requirements; `val::Tuple`, 
if the elements of the tuple fulfills those requirements.

Specify `domains` to only integrate over a part of the grid.
It should contain the names of the `AssemblyDomains` to integrate over.
Single, `String`, or collections; `Set{String}`, `AbstractVector{String}`, or `NTuple{N,String}`
inputs are supported. If `domains=nothing`, all domains are integrated. 

**Example:** *Calculate the average value and gradient* \\
We assume that we have done the setup: `buffer, state = setup_assembly(...)`,\\
assembled `K` and `r`, and solved `a=K\\r`
```julia
integrator = SimpleIntegrator((u, ∇u, state)->(1.0, u, ∇u), (0.0, 0.0, zero(Vec{dim})))
doassemble!(integrator, states, buffer; a=a)
area = integrator.val[1]
u_avg = integrator.val[2]/area 
∇u_avg = integrator.val[3]/area
```
If the `buffer` is setup to be threaded, this calculation will also be threaded. 
"""
mutable struct SimpleIntegrator{F<:Function, T, D<:Union{Nothing,Set{String}}}
    const fun::F
    val::T
    domains::D
end
SimpleIntegrator(fun::Function, val; domains=nothing) = SimpleIntegrator(fun, val, domains)
function SimpleIntegrator(fun::Function, val, domains::Union{AbstractVector{String},NTuple{<:Any,String}})
    return SimpleIntegrator(fun, val, Set(domains))
end
SimpleIntegrator(fun::Function, val, domain::String) = SimpleIntegrator(fun, val, Set((domain,)))

# Properties:
can_thread(::SimpleIntegrator) = true
need_colors(::SimpleIntegrator) = false
skip_this_domain(::SimpleIntegrator{<:Function,<:Any,Nothing}, ::String) = false
skip_this_domain(ig::SimpleIntegrator{<:Function,<:Any,<:Set}, name::String) = name ∉ ig.domains

function work_single_cell!(assembler::SimpleIntegrator, cellbuffer)
    cv = get_values(cellbuffer)
    cell_state = get_new_state(cellbuffer)
    ae = get_ae(cellbuffer)
    simple_integrate_cell!(assembler, cell_state, ae, cv, cellbuffer)
end

# Helper function to get state for qp for both vector and non-vector cell_state 
_get_qp_state(cell_state::AbstractVector, i) = cell_state[i]
_get_qp_state(cell_state, _) = cell_state # If not vector, just send the full. Up to user to know. 

# Updating of the value 
function add_to_val!(integrator::SimpleIntegrator{<:Function,<:Tuple}, f_val::Tuple, dΩ)
    integrator.val = map((v,f)-> v + f*dΩ, integrator.val, f_val)
end
function add_to_val!(integrator::SimpleIntegrator, f_val, dΩ)
    integrator.val += f_val*dΩ
end
zero_val!(integrator::SimpleIntegrator{<:Function,<:Tuple}) = (integrator.val = map(zero, integrator.val))
zero_val!(integrator::SimpleIntegrator) = (integrator.val = zero(integrator.val))

function simple_integrate_cell!(integrator::SimpleIntegrator, cell_state, ae, cv::CellValues, cellbuffer)
    length(Ferrite.getfieldnames(cellbuffer)) == 1 || throw(DimensionMismatch("Only one field supported for `CellValues`"))
    for q_point in 1:getnquadpoints(cv)
        dΩ = getdetJdV(cv, q_point)
        u = function_value(cv, q_point, ae)
        ∇u = function_gradient(cv, q_point, ae)
        f_val = integrator.fun(u, ∇u, _get_qp_state(cell_state, q_point))
        add_to_val!(integrator, f_val, dΩ)
    end
end

function simple_integrate_cell!(integrator::SimpleIntegrator, cell_state, ae, cv::NamedTuple, cellbuffer)
    length(Ferrite.getfieldnames(cellbuffer)) == length(cv) || throw(DimensionMismatch("Number of fields must match length of cellvalues tuple"))
    cv0 = first(values(cv))
    for q_point in 1:getnquadpoints(cv0)
        dΩ = getdetJdV(cv0, q_point)
        # Haven't benchmarked and `NamedTuple`s are tricky: Might need optimization to avoid dynamical dispatch/allocs:
        u = NamedTuple{keys(cv)}(map((k,v) -> function_value(v, q_point, ae, dof_range(cellbuffer, k)), keys(cv), values(cv)))
        ∇u = NamedTuple{keys(cv)}(map((k,v) -> function_gradient(v, q_point, ae, dof_range(cellbuffer, k)), keys(cv), values(cv)))
        f_val = integrator.fun(u, ∇u, _get_qp_state(cell_state, q_point))
        add_to_val!(integrator, f_val, dΩ)
    end
end

create_local(integrator::SimpleIntegrator) = SimpleIntegrator(integrator.fun, zero(integrator.val), integrator.domains)
create_local(integrator::SimpleIntegrator{<:Function,<:Tuple}) = SimpleIntegrator(integrator.fun, map(zero, integrator.val), integrator.domains)
scatter!(::SimpleIntegrator, ::SimpleIntegrator) = nothing
function gather!(base::TI, task::TI) where {TI<:SimpleIntegrator}
    add_to_val!(base, task.val, 1)
    zero_val!(task)
    return nothing
end
