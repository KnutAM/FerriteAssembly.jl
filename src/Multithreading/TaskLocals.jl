"""
    TaskLocals(base; num_tasks=Threads.nthreads())

TaskLocals is used to store a `base` value and a vector, `locals`. 
Items in `locals` are for each task and can be obtained with
[`get_local`](@ref). `base` can be obtained with [`get_base`](@ref).

A `base`-type that supports the `TaskLocals` API should define the following functions
- `create_local`
- `scatter!(local, base)`
- `gather!(base, local)`
"""
struct TaskLocals{TB,TL}
    base::TB
    locals::Vector{TL}
end
function TaskLocals(base; num_tasks=Threads.nthreads())
    locals = [create_local(base) for _ in 1:num_tasks]
    return TaskLocals(base, locals)
end

"""
    get_local(tl::TaskLocals, i::Int)

Get the `i`th `local` variable from `tl`
"""
get_local(tl::TaskLocals, i::Int) = tl.locals[i]

"""
    get_locals(tl::TaskLocals)

Get the vector of locals in `tl`
"""
get_locals(tl) = tl.locals

"""
    get_base(tl::TaskLocals)
    get_base(x::Any)

Get the `base` variable from `tl`, or, 
for an `x` that is not a `TaskLocals`,
return itself
"""
get_base(tl::TaskLocals) = tl.base
get_base(x) = x

function scatter!(tl::TaskLocals)
    map(Base.Fix2(scatter!, tl.base), get_locals(tl))
end

function gather!(tl::TaskLocals)
    map(Base.Fix1(gather!, tl.base), get_locals(tl))
end

# Each base::TB must implement the following methods
"""
    create_local(base::TB)

Create a task local variable, `local::TL`. In many cases, `TL=TB`.
`local`'s state should match that after calling `gather!(base, local)` 
(reset `local`) followed by `scatter!(local, base)` (add info from `base`)
(This is important, because it allows `create_local` to be called while working with 
tasks, if more locals are needed)
"""
function create_local end

"""
    scatter!(local, base)

Write any information from base required to be forwarded to 
the `local::typeof(create_local(base))`.
(typically if something in base have changed since its creation)
"""
function scatter! end

"""
    gather!(base, local)

Take any information from `local::typeof(create_local(base))` that should 
be added to `base` after an assembly. 
In addition, any accumulative information in `local` should 
be reset, such that running `scatter!` -> "do work" -> `gather!` should 
only affect the values in `base`. 
"""
function gather! end


# Idea for generic test interfacet
#=
function test_task_locals(modify_value!, read_value, base)
    # Test 
    # 1) modify_value! and read_value function supplied for testing works 
    # 2) That modifying base - creating task - gathering task does not change base
    #    as it is assumed that newly created task types have zero value 
    # 3) That modifying task will affect base after gather!

    vb0 = read_value(base)  # Get the initial value from base
    modify_value!(base)     # Modify the value in base
    vb1 = read_value(base)
    @test vb0 != vb1        # Check that value is updated (1)
    task = create_local(base)
    vt0 = read_value(task)  
    gather!(base, task)     
    vb2 = read_value(base)
    @test vb1 == vb2        # Check that gathering a newly created task doesn't contribute (2)
    modify_value!(task)
    vt1 = read_value(task)
    @test vt0 != vt1        # Check that value is updated (1)
    gather!(base, task) 
    @test vb2 != read_value(base) # Check that base was updated due to new value  (3)
end
=#
    
