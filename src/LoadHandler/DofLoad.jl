"""
    DofLoad(dofs::Union{Int, Vector{Int}}, load_function::Function)

Adds `load_function(t)` to `f[d]` for each degree of freedom number `d` in `dofs`.
Here, `t` is the simulation time passed to the [`LoadHandler`](@ref).
"""
struct DofLoad{F}
    dofs::Vector{Int}
    load_function::F # a[dofs] = load_function(t)
end
DofLoad(dof::Int, load_function) = DofLoad([dof], load_function)

function apply_dofload!(f::Vector, dof_load::DofLoad, time)
    val = dof_load.load_function(time)
    for dof in dof_load.dofs
        f[dof] += val
    end
end
