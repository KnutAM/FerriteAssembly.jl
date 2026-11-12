```@meta
CurrentModule = FerriteAssembly
```
# Algebraic variables

!!! warning "Experimental"
    The support for algebraic variables is experimental. The functions described on this page
    are not yet part of the public API, and may change without being considered a breaking change.

`Ferrite`'s algebraic variables (`AlgebraicVariable`) are global unknowns that are not associated
with the grid, for example a Lagrange multiplier that enforces the mean value of a field, a
cavity pressure, or the load factor in a path-following method. Their degrees of freedom are not
part of `celldofs`. To use them in the element or facet routines, they are given to the
[`DomainSpec`](@ref) for the domains that use them,
```julia
dh = DofHandler(grid)
add!(dh, :u, ip)
add!(dh, :λ, AlgebraicVariable())
close!(dh)
domain = DomainSpec(dh, material, cellvalues; algebraic_per_cell = (:λ,))
```
The same keyword is used for facet domains, and the same variable can be used in several domains,
including both cell and facet domains. The contributions from the different domains are then added
during assembly. See [Mean value constraints](@ref) for a complete example.

## Local degrees of freedom
For each item (cell or facet) in a domain with algebraic variables, the degrees of freedom of
each variable in `algebraic_per_cell` are added after the cell's degrees of freedom, in the given order.
Consequently, the following have the size of, and are ordered according to, these local degrees of freedom:
* `ae` and `get_aeold(buffer)`, including the current and old values of the algebraic variables
* `re` and `Ke`, including the contributions to the equations for the algebraic variables
* `celldofs(buffer)`, including the global degrees of freedom of the algebraic variables

The local dof range for each field is not changed, and `dof_range(buffer, name)` also gives the
range for an algebraic variable with name `name`. Since `ae` contains more values than the number
of shape functions, the element routine must pass the dof range when evaluating a field, e.g.
`function_value(cellvalues, q_point, ae, dof_range(buffer, :u))`.
The value of an algebraic variable is obtained as
```julia
av = FerriteAssembly.get_algebraic_values(buffer, :λ)
λ = algebraic_value(av, ae, dof_range(buffer, :λ))
```
which, e.g., gives a tensor for an `AlgebraicVariable{SymmetricTensor{2, 3}}`.
`Ferrite.getfieldnames(buffer)` and the `dofrange` given to [`create_cell_state`](@ref) contain the fields,
followed by the algebraic variables.

```@docs
FerriteAssembly.get_algebraic_values
```

## Sparsity pattern
The entries for the algebraic variables are not included in the default sparsity pattern. The required
coupling information for the domains is created by
```@docs
FerriteAssembly.algebraic_couplings
```
`Ferrite`'s `allocate_matrix` also accounts for constraints, e.g. affine constraints where the degrees of
freedom of an algebraic variable are the masters. Constraining the degrees of freedom of algebraic variables
is done with `AffineConstraint`s, and these can be applied both locally (with a `ConstraintHandler` given
to the [`KeReAssembler`](@ref)) and to the global system after assembly.
A degree of freedom that doesn't get any contributions, e.g. a component of a tensorial variable that
isn't used in any element, must be constrained to avoid a singular system matrix.

## Multithreading
The degrees of freedom of an algebraic variable are shared between all items in the domain. The coloring
used for multithreading does therefore not prevent race conditions when assembling. Consequently, a threaded
domain with algebraic variables is worked sequentially, unless the worker supports multithreading for
such domains:
```@docs
FerriteAssembly.can_thread_shared_dofs
```
For assembling the system matrix and residual, this is supported by `KeReAssembler(K, r; atomic = true)`,
where the task-local assemblers are atomic `Ferrite` assemblers. The result of such an assembly is not
deterministic to machine precision, since the order of the additions depends on the task scheduling.
The [`ReAssembler`](@ref), `Ferrite`'s assemblers, and a `KeReAssembler` created without `atomic = true`
work such domains sequentially. The [`SimpleIntegrator`](@ref) and [`QuadPointEvaluator`](@ref) work
multithreaded as for other domains.

## Further details
* **Automatic differentiation**: Defining only `element_residual!` works as for other cell domains,
  and the stiffness then includes the coupling to the algebraic variables.
  There is no automatic differentiation for facet domains.
* **Renumbering**: The global degrees of freedom of the algebraic variables are updated once
  when starting to work a domain, such that renumbering the dofhandler after
  setting up the domain buffer is supported.
* **Residual scaling**: [`ElementResidualScaling`](@ref) calculates one factor for each algebraic variable,
  in addition to the factor for each field.
* **Workers evaluating fields**: The [`SimpleIntegrator`](@ref) and [`QuadPointEvaluator`](@ref) pass the
  values and gradients of the fields. The values of the algebraic variables are not passed.
* **Loads**: The contributions added by the [`LoadHandler`](@ref) cannot depend on algebraic variables.
  Such contributions can be implemented with `facet_routine!` or `element_routine!`.
* **Coupled simulations**: A coupled buffer, see [`get_coupled_buffer`](@ref), gives access to the
  algebraic variables of the coupled simulation in the same way as for the buffer itself.
