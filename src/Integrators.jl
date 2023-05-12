# Integrators are a type of Assembler, but where the goal it to obtain one (or a few) values resulting from
# integrating a function over a domain. 

"""
    Integrator(fun::Function, val)

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

**Example:** *Calculate the average value and gradient* \\
We assume that we have done the setup: `buffer, state = setup_assembly(...)`,\\
assembled `K` and `r`, and solved `a=K\\r`
```julia
integrator = Integrator((u, ∇u, state)->(1.0, u, ∇u), (0.0, 0.0, zero(Vec{dim})))
doassemble!(integrator, states, buffer; a=a)
area = integrator.val[1]
u_avg = integrator.val[2]/area 
∇u_avg = integrator.val[3]/area
```
If the `buffer` is setup to be threaded, this calculation will also be threaded. 
"""
mutable struct Integrator{F<:Function, T}
    const fun::F
    val::T
end
function assemble_cell_reinited!(assembler::Integrator, cell_state, cellbuffer)
    cv = get_cellvalues(cellbuffer)
    ae = get_ae(cellbuffer)
    integrate_cell!(assembler, cell_state, ae, cv, cellbuffer)
end

# Helper function to get state for qp for both vector and non-vector cell_state 
_get_qp_state(cell_state::AbstractVector, i) = cell_state[i]
_get_qp_state(cell_state, _) = cell_state # If not vector, just send the full. Up to user to know. 

# Updating of the value 
function update_val!(integrator::Integrator{<:Function,<:Tuple}, f_val::Tuple, dΩ)
    integrator.val = map((v,f)-> v + f*dΩ, integrator.val, f_val)
end
function update_val!(integrator::Integrator, f_val, dΩ)
    integrator.val += f_val*dΩ
end

function integrate_cell!(integrator::Integrator, cell_state, ae, cv::CellValues, cellbuffer)
    length(Ferrite.getfieldnames(cellbuffer)) == 1 || throw(DimensionMismatch("Only one field supported for `CellValues`"))
    for q_point in 1:getnquadpoints(cv)
        dΩ = getdetJdV(cv, q_point)
        u = function_value(cv, q_point, ae)
        ∇u = function_gradient(cv, q_point, ae)
        f_val = integrator.fun(u, ∇u, _get_qp_state(cell_state, q_point))
        update_val!(integrator, f_val, dΩ)
    end
end

function integrate_cell!(integrator::Integrator, cell_state, ae, cv::NamedTuple, cellbuffer)
    length(Ferrite.getfieldnames(cellbuffer)) == length(cv) || throw(DimensionMismatch("Number of fields must match length of cellvalues tuple"))
    cv0 = first(values(cv))
    for q_point in 1:getnquadpoints(cv0)
        dΩ = getdetJdV(cv0, q_point)
        # Haven't benchmarked and `NamedTuple`s are tricky: Might need optimization to avoid dynamical dispatch/allocs:
        u = NamedTuple{keys(cv)}(map((k,v) -> function_value(v, q_point, ae, dof_range(cellbuffer, k)), keys(cv), values(cv)))
        ∇u = NamedTuple{keys(cv)}(map((k,v) -> function_gradient(v, q_point, ae, dof_range(cellbuffer, k)), keys(cv), values(cv)))
        f_val = integrator.fun(u, ∇u, _get_qp_state(cell_state, q_point))
        update_val!(integrator, f_val, dΩ)
    end
end

create_local(integrator::Integrator) = Integrator(integrator.fun, zero(integrator.val))
create_local(integrator::Integrator{<:Function,<:Tuple}) = Integrator(integrator.fun, map(zero, integrator.val))
scatter!(::Integrator, ::Integrator) = nothing
function gather!(base::TI, task::TI) where {TI<:Integrator}
    update_val!(base, task.val, 1)
    update_val!(task, task.val, 0) # Makes it work also for Tuple, and speed is not critical as not called so often. 
    return nothing
end