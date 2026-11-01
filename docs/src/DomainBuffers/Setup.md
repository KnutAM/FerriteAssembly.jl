```@meta
CurrentModule = FerriteAssembly
```
# [Setup](@id DomainBufferSetup)
## Setup API
```@docs
DomainSpec
setup_domainbuffer
setup_domainbuffers
```

## AbstractDomainBuffer
The domain buffer be a `DomainBuffer`, `ThreadedDomainBuffer`, or a [`FerriteAssembly.DomainBuffers`](@ref)
(an immutable `AbstractDict{String}`) with eltype of one of the former. The
following functions are defined for these buffers:
```@docs
FerriteAssembly.DomainBuffers
FerriteAssembly.get_material(::FerriteAssembly.DomainBuffers, ::String)
FerriteAssembly.get_dofhandler(::FerriteAssembly.DomainBuffers)
FerriteAssembly.get_grid
FerriteAssembly.get_state(::FerriteAssembly.DomainBuffers, ::String)
FerriteAssembly.get_old_state(::FerriteAssembly.DomainBuffers, ::String)
FerriteAssembly.getset
update_states!(::FerriteAssembly.DomainBuffers)
revert_states!(::FerriteAssembly.DomainBuffers)
set_time_increment!(::FerriteAssembly.DomainBuffers, ::Any)
```

## Updating materials
!!! warning "Material update contract"
    `get_material` (above) always returns the *base* material; a threaded domain buffer's
    task-local copies, used by threaded `work!`, are **not** kept in sync with mutations to
    that returned object. Use `replace_material` instead to create a new buffer if this is required.
```@docs
FerriteAssembly.replace_material(::FerriteAssembly.DomainBuffers, ::Any)
FerriteAssembly.replace_material(::FerriteAssembly.DomainBuffers, ::String, ::Any)
```

## Concurrent `work!` calls
Multiple `work!` calls — including threaded ones — may run concurrently (e.g. each spawned
with `Threads.@spawn`, as when running several points of a parameter sweep in parallel),
*as long as each call's `Simulation`/domain buffer owns genuinely independent mutable
storage*: its own item-buffer scratch (`Ke`, `re`, `ae`, ...), its own state variables, and
its own output (`K`, `r`, or whatever else the `worker` writes into). Two calls that share
any of that storage must not run concurrently with each other. In particular:
* Working the same buffer object (or `Simulation`) twice at once is not safe.
* A buffer and a sibling produced from it with `replace_material` (or vice versa) share
  storage: `replace_material` only swaps out the `material` field, so everything else
  (item-buffer scratch and state variables) remains the *same* object, by reference. Give
  each concurrently-worked variant (e.g. each parameter-sweep point) its own, independently
  `setup_domainbuffer`'d buffer instead of deriving it with `replace_material`.
* Members of a [`CoupledSimulations`](@ref) group that share storage must not be worked
  concurrently — see the "Concurrency contract" warning below.
* `user_data` (see `DomainSpec`) is passed by reference; sharing the same `user_data`
  between two buffers worked concurrently is only safe if neither `worker` mutates it.

When running many `work!` calls concurrently this way, each threaded buffer's own
`num_tasks` (defaulting to `Threads.nthreads()`) is set independently of how many *outer*
calls happen to run at once: `P` concurrent outer calls spawn up to `P * num_tasks` tasks
in total. This is not a correctness problem — Julia's scheduler tolerates oversubscription —
but it is wasted overhead. Set `num_tasks` explicitly per buffer (e.g.
`max(1, Threads.nthreads() ÷ P)`) to avoid spawning far more tasks than there are threads
available to run them.

## Coupled simulations
The `Simulation` type contains an abstract domain buffer, along with (optionally) 
the global degree of freedom values, which are used to get the local values for each item.

A [`CoupledSimulations`](@ref) group is built, once, from a set of named `Simulation`s, and
gives access to values from other simulations at the item level (e.g. state variables and
local dof-values) via [`get_coupled_buffer`](@ref). For example, when solving two separate
problems in parallel using staggered iterations. See the
[Phase-field fracture tutorial](@ref Phase-field-fracture) for an example. Coupling is
resolved entirely at group-construction time; `work!`ing a group member (`work!(worker,
group.member_name)`) never re-discovers or rebuilds the coupling.

!!! warning "Concurrency contract"
    Coupled buffers reference the partner's *actual* mutable storage — no copies are made.
    `work!` calls that share any of that storage must therefore not run concurrently with
    each other. This includes: working two members of the *same* group at the same time;
    working a member of a group at the same time as its own original (pre-group) source
    `Simulation`; working members of *two different* groups that were built from the same
    source `Simulation`(s) (e.g. a group and a later `replace_material`-built group that
    still shares some members' storage by reference); and re-entrant `work!` calls that
    would reuse the same scratch. Ordinary staggered iteration — working one member, then
    another, in sequence — is safe; it is *simultaneous* access to shared scratch that is not.
```@docs
Simulation
CoupledSimulations
FerriteAssembly.replace_material(::CoupledSimulations, ::Symbol, ::Any)
```
