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
domain = DomainSpec(dh, material, cellvalues; algebraic_variables = (:λ,))
```
The same keyword is used for facet domains, and the same variable can be used in several domains,
including both cell and facet domains. The contributions from the different domains are then added
during assembly. See [Mean value constraints](@ref) for a complete example.

Each element in `algebraic_variables` is either the name of a *standalone* variable, as above, or
an [`AttachAlgebraic`](@ref FerriteAssembly.AttachAlgebraic) declaration of a variable *attached*
to a field, see [Attached variables](@ref).

## Local degrees of freedom
For each item (cell or facet) in a domain with standalone algebraic variables, the degrees of freedom of
each standalone variable are added after the cell's degrees of freedom, in the given order.
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
followed by the algebraic variables. If not all components of a variable are used in a domain (see
[`AttachAlgebraic`](@ref FerriteAssembly.AttachAlgebraic)), its value cannot be obtained this way.

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
Multithreaded assembly of domains with algebraic variables requires an atomic assembler, e.g.
`KeReAssembler(K, r; atomic = true)`. Otherwise, such domains are assembled sequentially.

## Attached variables
An algebraic variable can also be *attached* to a field. The degrees of freedom of the variable then act
as additional degrees of freedom of the field, which is useful when the field contains a global contribution,
for example
* the constant out-of-plane strain, ``\varepsilon_{33}``, in generalized plane strain, see
  [Generalized plane stress](@ref);
* the macroscopic displacement gradient, ``\boldsymbol{H}``, in computational homogenization with
  ``\boldsymbol{u} = \boldsymbol{H}\cdot[\boldsymbol{x} - \boldsymbol{x}_\mathrm{c}] + \tilde{\boldsymbol{u}}``;
* rigid body modes or global enrichment modes.

The attachment is declared with
```@docs
FerriteAssembly.AttachAlgebraic
```
for example,
```julia
dh = DofHandler(grid)
add!(dh, :u, ip^2)
add!(dh, :p, ip)
add!(dh, :H, AlgebraicVariable{Tensor{2, 2}}())
add!(dh, :λ, AlgebraicVariable())
close!(dh)
domain = DomainSpec(dh, material, (u = cellvalues_u, p = cellvalues_p);
    algebraic_variables = (FerriteAssembly.AttachAlgebraic(:H => :u; indices = ((1, 2), (2, 1))), :λ))
```
For each field, in the order of the `SubDofHandler`, the local degrees of freedom are the field's degrees of
freedom, followed by those of the variables attached to the field, in declaration order. Then the standalone
variables follow. In the example above, with linear quadrilaterals, the local dof ranges are
`u = 1:10`, `H = 9:10`, `p = 11:14`, and `λ = 15:15`. Hence, `dof_range(buffer, :u)` includes the attached
degrees of freedom, such that the values for `:u` (`cellvalues_u`) must provide 10 base functions, where
the last two correspond to the attached components of `H`. Without standalone variables, an element routine
written for a problem without algebraic variables can therefore be used unchanged, as long as the values
for the attached field provide the additional base functions. This is checked when setting up the domain.
For the field's own degrees of freedom, use
```@docs
FerriteAssembly.field_dof_range
```

Note that the [`LoadHandler`](@ref) only adds contributions to the field's own degrees of freedom. If an
external load does work on the attached base functions, this must be included in the element routine, or added to
the global residual at `algebraic_dofs(dh, name)`.

### Values for attached variables
`FerriteAssembly` does not provide values types with additional base functions; these are problem specific.
A custom values type for a cell domain should subtype `Ferrite.AbstractCellValues` and define
* `getnbasefunctions`, including the attached base functions;
* `getnquadpoints` and `getdetJdV`;
* `shape_value` and `shape_gradient` for all base functions. These must return the same type for all
  base functions, and should not allocate. The value and the gradient are independent, i.e. the gradient of
  an attached base function need not be the gradient of its value (e.g. a zero value with a non-zero
  out-of-plane gradient);
* `getngeobasefunctions` and `geometric_value` (or `spatial_coordinate`);
* `reinit!(values, cell::Union{Nothing, Ferrite.AbstractCell}, x::AbstractVector{<:Vec})`, where data
  depending on the cell can be updated (e.g. the position of the quadrature points). `deepcopy` must give
  an independent copy, as used for multithreading;
* optionally `shape_hessian`, if used by the element routine.
`Ferrite`'s generic methods then give the remaining accessors, e.g. `shape_symmetric_gradient`,
`shape_divergence`, `function_value`, and `function_gradient`. For a facet domain, the values type should
subtype `Ferrite.AbstractFacetValues` and define `reinit!(values, cell, x, facetnr)` instead.
A custom values type can be checked with
```@docs
FerriteAssembly.check_values_interface
```
If the values type lifts a problem to a higher dimension, e.g. to use a 3d element routine for a 2d grid,
the coordinates given to the element routine can be lifted as well, by overloading
```@docs
FerriteAssembly.routine_coordinates
```

### Attached and standalone variables in the same domain
The degrees of freedom of the standalone variables follow those of the fields, such that `ae`, `re`, and `Ke`
contain more entries than the base functions of the fields. To still use an existing routine for the fields,
including their attached variables, it can be called via
```@docs
FerriteAssembly.field_block
```

## Further details
* **Automatic differentiation**: Defining only `element_residual!` works as for other cell domains,
  and the stiffness then includes the coupling to the algebraic variables.
  There is no automatic differentiation for facet domains.
* **Renumbering**: The global degrees of freedom of the algebraic variables are updated once
  when starting to work a domain, such that renumbering the dofhandler after
  setting up the domain buffer is supported.
* **Residual scaling**: [`ElementResidualScaling`](@ref) calculates one factor for each algebraic variable,
  in addition to the factor for each field. The factor for a field only includes the field's own degrees
  of freedom, see [`FerriteAssembly.field_dof_range`](@ref).
* **Workers evaluating fields**: The [`SimpleIntegrator`](@ref) and [`QuadPointEvaluator`](@ref) pass the
  values and gradients of the fields, including the contributions from attached variables. The values of
  the standalone algebraic variables are not passed.
* **Loads**: The contributions added by the [`LoadHandler`](@ref) cannot depend on algebraic variables.
  Such contributions can be implemented with `facet_routine!` or `element_routine!`.
* **Coupled simulations**: A coupled buffer, see [`get_coupled_buffer`](@ref), gives access to the
  algebraic variables of the coupled simulation in the same way as for the buffer itself.
