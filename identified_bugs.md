# Identified bugs and inconsistencies

Audit date: 2026-09-08. Reviewed commit: `5f56bfdae413631cd1d3817648f2ac8fb97a9c4e` (package version 0.3.7).

This is a work queue, ordered by severity. All issues are open; check an item off after implementing its correction and regression check. IDs should remain stable when issues are resolved. Locations refer to the audited revision.

## Scope and verification

Reviewed all package source files, all test files, the tutorials and how-tos, reference documentation, documentation tooling, project/manifest files, and CI workflows. Existing untracked `CLAUDE.md` and `debug.md` were left unchanged. Diagnostic claims in `debug.md` were checked against current implementation rather than treated as findings by themselves. No implementation fixes were made.

The complete `test/runtests.jl` suite passed **1,689 assertions** with Julia 1.12.5 and four worker threads, using the existing documentation environment (Ferrite 1.1.0, MaterialModelsBase 0.4.0). Command, with the installed Julia binary substituted for the launcher:

```sh
julia --compiled-modules=existing --project=docs --threads=4 -e 'include("test/runtests.jl")'
```

Additional isolated probes reproduced the failures labeled **Reproduced** below. The threaded-assembly and local-constraints how-tos were also executed. The affine-constraint stress probe disagreed with sequential assembly on all 20 repetitions.

Limitations: this was not a fresh dependency resolution or a run of every supported dependency/OS combination. Julia 1.11.9 could not load the existing Julia-1.12-resolved docs environment (`PrecompileTools` referenced `Base.StaticData`); that environment failure is not evidence that the package itself fails on a freshly resolved Julia 1.11 environment. The full documentation build, mesh downloads, and all complete tutorial simulations were not executed. Findings in those paths are labeled with the narrower evidence actually obtained. An audit cannot establish the absence of further bugs.

## Work queue

| Done | ID | Severity | Issue |
| --- | --- | --- | --- |
| [x] | BUG-001 | High | Affine constraints race during threaded assembly |
| [x] | BUG-002 | High | Empty custom chunks silently terminate workers |
| [ ] | BUG-003 | High | Coupled threaded buffers read stale time increments |
| [ ] | BUG-004 | High | Mixed-material tutorial accumulates external loads between steps |
| [ ] | BUG-005 | High | Mixed-material tutorial computes incorrect plane-stress output |
| [ ] | BUG-006 | Medium | Mutating the public material affects sequential and threaded work differently |
| [ ] | BUG-007 | Medium | Nonpositive task counts silently disable all work |
| [ ] | BUG-008 | Medium | Threaded load handlers fail when either distributed-load category is empty |
| [ ] | BUG-009 | Medium | Accepted abstract solution vectors fail during buffer initialization |
| [ ] | BUG-010 | Medium | Public domain-dictionary type is not accepted by work dispatch |
| [ ] | BUG-011 | Medium | Mixed cell/facet setup depends on dictionary iteration order |
| [ ] | BUG-012 | Medium | Documented quadrature-rule load inputs have no implementation |
| [ ] | BUG-013 | Medium | Concrete Grid dispatch prevents threaded IGA/custom-grid use |
| [ ] | BUG-014 | Medium | Function-based quadrature evaluation rejects whole-cell states |
| [ ] | BUG-015 | Medium | Threaded-assembly how-to produces a NaN residual |
| [ ] | BUG-016 | Low | Empty cell domains fail with internal indexing errors |
| [ ] | BUG-017 | Low | Coupling setup incompletely validates domain/task compatibility |
| [ ] | BUG-018 | Low | Public state accessors promise dictionaries but return a limited wrapper |
| [ ] | BUG-019 | Low | Assembler documentation promises nonexistent facet autodiff |
| [ ] | BUG-020 | Low | WeakForm examples use obsolete loading APIs and the wrong field/sign |
| [ ] | BUG-021 | Low | Local-constraint example compares a solution against itself |
| [ ] | BUG-022 | Low | Some tests labeled as threaded actually run sequentially |
| [ ] | BUG-023 | Low | Figure-factory manifest is incompatible with current package source |

For the small reproductions below, this common setup creates two scalar quadrilateral cells with total area 4:

```julia
using Ferrite, FerriteAssembly
import FerriteAssembly as FA

grid = generate_grid(Quadrilateral, (2, 1))
ip = Lagrange{RefQuadrilateral, 1}()
dh = close!(add!(DofHandler(grid), :u, ip))
cv = CellValues(QuadratureRule{RefQuadrilateral}(2), ip)
db = setup_domainbuffer(DomainSpec(dh, nothing, cv))
```

## BUG-001 — Affine constraints race during threaded assembly

**Severity:** High. **Evidence:** Reproduced, including numerical disagreement on 20/20 runs.

**Location:** `src/Workers/Assemblers.jl:139–154` (`can_thread`, `assemble_contributions!`); `src/work.jl:64–89`.

**Problem:** `KeReAssembler` always allows threading, including when `ch` contains affine constraints. Its local assembly calls `Ferrite.apply_assemble!`, which can write to master degrees of freedom outside the current cell. Ordinary mesh coloring does not protect these extra writes. Task-local assembler scratch does not solve this: the race is on the shared global matrix/vector.

**Failure scenario:** On a 2,000-cell line grid, constrain every degree of freedom except DOF 1 to DOF 1 with coefficient 1. Use an element routine that fills its local matrix and residual with ones. Allocate the global matrix with `allocate_matrix(dh, ch)`. Sequential assembly gives `(K[1,1], r[1]) == (8000.0, 4000.0)`. Four-task assembly produced values such as `(5897.0, 3320.0)` and `(4989.0, 3004.0)`, varying between runs, with no error.

**Smallest reasonable correction:** Make `can_thread` return false when the constraint handler contains nontrivial affine constraints, using the existing sequential fallback. Preserve threading for cases proven safe, such as ordinary cell-local Dirichlet application. Constraint-aware coloring or synchronized assembly is a larger alternative.

**Regression check:** Add an affine-constraint case whose otherwise independent cells share master DOFs. Compare both matrix and residual with sequential assembly, and assert the chosen fallback for unsafe constraints. Existing `test/assemblers.jl` tests only sequential local constraints.

## BUG-002 — Empty custom chunks silently terminate workers

**Severity:** High. **Evidence:** Reproduced deterministically with one task.

**Location:** `src/Multithreading/TaskChunks.jl:25–32,83–99`; `src/work.jl:77–78`.

**Problem:** An empty vector represents both a legitimate empty chunk and the end of the queue. Custom chunk validation allows empty chunks anywhere. A worker exits when it retrieves one, even if later chunks contain work. With enough leading empty chunks, all workers exit.

**Reproduction:** Using the common setup:

```julia
b = setup_domainbuffer(
    DomainSpec(dh, nothing, cv; chunks = [[Int[], [1, 2]]]);
    threading = true, num_tasks = 1)
ig = SimpleIntegrator(Returns(1.0), 0.0)
work!(ig, b)
ig.val # 0.0; expected 4.0
```

**Smallest reasonable correction:** Remove/reject empty chunks at setup, or distinguish exhaustion with a separate sentinel and skip empty chunks. Also avoid generating unnecessary empty chunks in `split_in_chunks`.

**Regression check:** Exercise empty chunks at the beginning, middle, and end, with one and several tasks. Assert that every cell is visited exactly once, not merely that the union of supplied chunk indices is correct.

## BUG-003 — Coupled threaded buffers read stale time increments

**Severity:** High. **Evidence:** Reproduced; sequential values were finite while threaded values were NaN.

**Location:** `src/DomainBuffers.jl:232–234,251–256`; `src/ItemBuffers/AbstractItemBuffer.jl:108–110`; `src/ItemBuffers/CellBuffer.jl:129–153`; `src/work.jl:64–67`.

**Problem:** `set_time_increment!` updates only a domain's base buffer. At work time, only the active simulation's task-local buffers are scattered. The buffers linked to coupled simulations are reinitialized for cell values and states, but their time increment is never refreshed. They retain the initial NaN or the value from the last time that partner was worked independently. Coupling also copies buffer objects, so relying on unrelated future scatters is insufficient in general.

**Failure scenario:** Create threaded domains A and B, link A to B, set B's increment to 0.25, then work A with `CoupledSimulations(b = Simulation(B))`. An element reading `get_time_increment(get_coupled_buffer(buffer, :b))` sees NaN. The same sequence with sequential domains sees 0.25.

**Smallest reasonable correction:** Refresh coupled-buffer time increments from the corresponding supplied simulation during coupled reinitialization, or explicitly scatter all linked coupled buffers before work. The source must be the current simulation's base value, not a potentially stale linked copy.

**Regression check:** Read the partner's increment before its first assembly and after changing it between staggered iterations. Check sequential/threaded and ordinary/autodiff buffers, including independently coupled copies as in the fracture tutorial.

## BUG-004 — Mixed-material tutorial accumulates external loads between steps

**Severity:** High. **Evidence:** Static control-flow evidence plus a reproduced load-history calculation.

**Location:** `docs/src/literate_tutorials/mixed_materials.jl:106,110–120`; additive behavior in `src/LoadHandler/LoadHandler.jl:42–55`.

**Problem:** `fext` is zeroed once, outside the time loop. `apply!(fext, lh, t)` adds the current load without clearing previous loads. The solver therefore applies the sum of all previous load levels instead of the load prescribed at the current time.

**Failure scenario:** The supplied 20-point history `range(0, 1, 20)` with traction proportional to time ends with ten times the intended final load. A small Neumann-load reproduction gave exactly `10.0` for the ratio of the accumulated final load to the correct final load. Plastic deformation and the published stress/displacement results consequently correspond to the wrong loading path.

**Smallest reasonable correction:** Call `fill!(fext, 0)` immediately before each time-step `apply!`. Regenerate numerical reference results and images affected by the loading change.

**Regression check:** Compare each time-step external vector with a newly allocated vector evaluated at that time. Include a repeated load level and unloading, where unintended accumulation is especially clear.

## BUG-005 — Mixed-material tutorial computes incorrect plane-stress output

**Severity:** High. **Evidence:** Reproduced with the tutorial's elastic material and postprocessing formula.

**Location:** `docs/src/literate_tutorials/mixed_materials.jl:53–60,84–90` (`calculate_stress`).

**Problem:** The materials are configured for `PlaneStress`, but postprocessing expands the in-plane strain with zero out-of-plane components and applies the 3D elastic tensor directly. Plane stress requires the eliminated out-of-plane strains to satisfy the out-of-plane stress constraint; setting those strains to zero instead computes a plane-strain-like response. Simply reducing the resulting stress tensor cannot correct it. The plastic path uses the same incorrect strain expansion.

**Failure scenario:** For the tutorial's `E = 210e3`, `ν = 0.3` elastic material, with `ε11 = 0.01`, `ε22 = ε12 = 0`, the tutorial formula gives `(σ11, σ22) = (2826.9231, 1211.5385)`. `MaterialModelsBase.material_response` for the configured plane-stress material gives `(2307.6923, 692.3077)`. Thus a successful solve can still export substantially incorrect stresses.

**Smallest reasonable correction:** Obtain stresses through a constitutively consistent reduced-stress evaluation, or retain the converged stress/full strain during assembly for postprocessing. For plastic materials, preserve the correct old/converged state semantics rather than inadvertently advancing the material again while plotting.

**Regression check:** Compare postprocessed stresses with those returned during assembly for elastic and plastic plane-stress points. Include an elastic analytical reference and a check on out-of-plane stress. Correct BUG-004 before regenerating tutorial results.

## BUG-006 — Mutating the public material affects sequential and threaded work differently

**Severity:** Medium. **Evidence:** Reproduced.

**Location:** `src/ItemBuffers/CellBuffer.jl:74–77`; `src/ItemBuffers/FacetBuffer.jl:86–90`; `src/ItemBuffers/AbstractItemBuffer.jl:108–110`; `src/DomainBuffers.jl:225`.

**Problem:** Task-local materials are deep-copied at setup. The public `get_material(db)` returns the base material, while `scatter!` forwards only the time increment. Changes to a mutable material or a mutable payload therefore affect sequential work, but leave threaded work using its original copies. The getter documentation does not explain this distinction.

**Failure scenario:** Set up a material with `k = 1.0`, then set `FA.get_material(db).k = 2.0`. A residual routine filling `re` with `m.k` uses 2.0 sequentially and 1.0 in threaded work. A sequential-only worker over the same threaded domain also observes the base material, creating inconsistent results within one simulation.

**Smallest reasonable correction:** Define and document an explicit material-update contract. Either synchronize supported material updates, or require `replace_material` and clearly document that mutating the getter's result does not propagate to task-local copies. Avoid blindly overwriting material/cache state that elements intentionally maintain per task.

**Regression check:** Mutate a parameter between two work calls and check the chosen supported update API in both modes; include a sequential-only worker operating on a threaded domain.

## BUG-007 — Nonpositive task counts silently disable all work

**Severity:** Medium. **Evidence:** Reproduced for `num_tasks = 0`; the same empty-range path handles negative counts.

**Location:** `src/DomainBuffers.jl:199–204`; `src/Multithreading/TaskLocals.jl:19–22`; `src/work.jl:72`.

**Problem:** No constructor checks that `num_tasks` is positive. Empty task-local arrays are created and the spawn loop executes zero times; `work!` returns successfully without visiting any item.

**Reproduction:** `setup_domainbuffer(DomainSpec(dh, nothing, cv); threading=true, num_tasks=0)` followed by a constant `SimpleIntegrator` returns 0.0 instead of area 4.0.

**Smallest reasonable correction:** Reject task counts below 1 with `ArgumentError` at setup, before creating buffers.

**Regression check:** Reject zero and negative counts and verify correct results for 1 and a count greater than the number of available threads.

## BUG-008 — Threaded load handlers fail when either distributed-load category is empty

**Severity:** Medium. **Evidence:** Reproduced.

**Location:** `src/LoadHandler/LoadHandler.jl:42–55`; `src/work.jl:31–33`; `src/DomainBuffers.jl:208`.

**Problem:** Applying a load handler always works both its Neumann and body-load dictionaries. Threaded multidomain work calculates `maximum(get_num_tasks, values(dbs))` even for an empty dictionary, which throws. Empty categories are normal for body-only, boundary-only, DOF-only, and entirely empty load handlers.

**Reproduction:** Using the common setup:

```julia
lh = LoadHandler(dh; threading=true)
add!(lh, DofLoad(1, t -> 2.0))
apply!(zeros(ndofs(dh)), lh, 1.0)
# ArgumentError: reducing over an empty collection is not allowed
```

**Smallest reasonable correction:** Make multidomain threaded `work!` return immediately for empty domains, before determining the task count. This also repairs callers other than `LoadHandler`.

**Regression check:** Run all four empty-category combinations above in both modes. The existing threaded load test always populates both dictionaries.

## BUG-009 — Accepted abstract solution vectors fail during buffer initialization

**Severity:** Medium. **Evidence:** Reproduced with a vector view.

**Location:** `src/Simulation.jl:10–17`; `src/Utils/utils.jl:31–38`; cell/facet calls to `_copydofs!`.

**Problem:** `Simulation` explicitly accepts `AbstractVector` for `a` and `aold`, but `_copydofs!` only accepts a concrete `Vector`. A validly constructed simulation therefore fails when work begins. Views are especially useful for passing components of a coupled global solution.

**Reproduction:** `work!(SimpleIntegrator((u,g,s)->u, 0.0), Simulation(db, view(ones(ndofs(dh)), :)))` throws `MethodError` for `_copydofs!`.

**Smallest reasonable correction:** Generalize the global input of `_copydofs!` to supported `AbstractVector` types, retaining bounds checks and DOF indexing semantics.

**Regression check:** Pass views as `a` and `aold` for cells and facets, sequentially and threaded. Compare with copies of the same views.

## BUG-010 — Public domain-dictionary type is not accepted by work dispatch

**Severity:** Medium. **Evidence:** Reproduced even for a dictionary containing one ordinary domain.

**Location:** `src/DomainBuffers.jl:3`; `src/Simulation.jl:21–24`; `src/work.jl:1–2,21–53`.

**Problem:** The public entry point accepts `Dict{String,<:AbstractDomainBuffer}`, but forwards to `Simulation` methods restricted to dictionaries whose declared value type is a subtype of `DomainBuffer` or `ThreadedDomainBuffer`. Widening the dictionary to the documented abstract type makes dispatch fail regardless of its actual contents. A dictionary mixing sequential and threaded buffers has the same problem.

**Reproduction:** `work!(SimpleIntegrator(Returns(1.0), 0.0), Dict{String,FA.AbstractDomainBuffer}("a" => db))` throws `MethodError` for `work!(..., ::Simulation{Dict{String,AbstractDomainBuffer},...})`.

**Smallest reasonable correction:** Add a generic multidomain simulation path that dispatches on each actual domain, preserving domain skipping and coupled-simulation routing. Keep the specialized homogeneous fast paths if beneficial.

**Regression check:** Test a widened dictionary with all-sequential values, all-threaded values, and mixed values. Each should agree with working the domains separately.

## BUG-011 — Mixed cell/facet setup depends on dictionary iteration order

**Severity:** Medium. **Evidence:** Reproduced in the cell-first validation branch.

**Location:** `src/setup.jl:61–89` (`check_input`).

**Problem:** Validation chooses its behavior solely from the first domain's item type. If it sees a cell domain first, it subsequently uses every domain's items as integer cell-array indices, including `FacetIndex` values. If it sees a facet domain first, it performs no checks. Thus otherwise usable cell and boundary contributions can fail during setup depending on dictionary keys/order. Cell counts also incorrectly include facet counts.

**Failure scenario:** Build one cell domain and one facet domain and call `FA.check_input(dbs, Int)`, the branch used when cells come first. It warns, then throws `ArgumentError: invalid index: FacetIndex(...)`. Suppressing setup warnings bypasses the failing validation rather than repairing it.

**Smallest reasonable correction:** Validate cell coverage only across cell domains, and validate other domain types separately. Do not choose the policy from one arbitrary dictionary entry.

**Regression check:** Set up mixed cell/facet dictionaries under several key arrangements and verify identical behavior and combined assembly results.

## BUG-012 — Documented quadrature-rule load inputs have no implementation

**Severity:** Medium. **Evidence:** Reproduced for both volume and facet rules.

**Location:** `src/LoadHandler/BodyLoad.jl:2–3,23–28`; `src/LoadHandler/Neumann.jl:2,27–32`; `src/LoadHandler/defaultvalues.jl:13–16,28–31`.

**Problem:** Load constructor documentation advertises quadrature-rule objects. The automatic values constructors only implement an integer order or an existing values object. Adding a load with a rule object fails. The Neumann documentation additionally calls the input a `QuadratureRule`, although boundary integration normally needs `FacetQuadratureRule`.

**Reproduction:** `add!(LoadHandler(dh), BodyLoad(:u, QuadratureRule{RefQuadrilateral}(2), Returns(1.0)))` throws `MethodError` for `autogenerate_cellvalues`. Passing `FacetQuadratureRule{RefQuadrilateral}(2)` to `Neumann` likewise fails for `autogenerate_facetvalues`.

**Smallest reasonable correction:** Implement construction from the appropriate volume/facet quadrature-rule objects and correct the signatures in the documentation, or remove the unsupported promise explicitly.

**Regression check:** Compare integer-order, explicit-rule, and explicit-values load vectors, including a rule that cannot be represented merely by selecting an integer order.

## BUG-013 — Concrete Grid dispatch prevents threaded IGA/custom-grid use

**Severity:** Medium. **Evidence:** Reproduced with the installed IGA package.

**Location:** `src/Multithreading/TaskChunks.jl:83,103,109,117`; `src/DomainBuffers.jl:199–203`; `docs/src/literate_tutorials/iga.jl`.

**Problem:** All `create_chunks` methods require a concrete `Ferrite.Grid`, while the surrounding DOF/buffer APIs support other grid implementations. Even the overload that only validates explicitly supplied chunks unnecessarily requires `Grid`. An IGA `BezierGrid` works sequentially but cannot enable threading, even with user-supplied chunks.

**Failure scenario:** Build a small plate-with-hole `BezierGrid`, an `IGAInterpolation`, and `BezierCellValues` as in the tutorial. Sequential `setup_domainbuffer` succeeds. Adding `threading=true` throws `MethodError: no method matching create_chunks(::BezierGrid, ::Vector{Int64}, ::Nothing)`.

**Smallest reasonable correction:** Generalize grid-independent chunk validation/conversion to the abstract grid interface. Use automatic coloring only where supported; otherwise allow valid supplied chunks and issue a clear limitation for automatic coloring.

**Regression check:** Exercise a non-`Grid` implementation with explicit safe chunks. Verify numerical agreement with sequential assembly; separately test whether its automatic coloring is supported.

## BUG-014 — Function-based quadrature evaluation rejects whole-cell states

**Severity:** Medium. **Evidence:** Reproduced with a scalar cell state.

**Location:** `src/Workers/QuadPointEvaluator.jl:86–99`; compare `src/Workers/Integrators.jl:125–127` (`_get_qp_state`).

**Problem:** The package supports both a whole-cell state and a vector of quadrature-point states. `SimpleIntegrator` handles both, but both function-based evaluator overloads require `cell_state::AbstractVector`. A callback that does not even use state still cannot evaluate a valid domain with whole-cell state. This restriction is absent from the evaluator's public documentation.

**Failure scenario:** Define `FA.create_cell_state(::MyMaterial, args...) = 1.0`, create the domain, then construct `QuadPointEvaluator{Float64}(db, (m,u,g,s)->s)`. Construction succeeds; `work!` fails with `MethodError` for `eval_quadpoints_cell!`.

**Smallest reasonable correction:** Generalize the state argument and reuse `_get_qp_state`, matching `SimpleIntegrator`. Alternatively, explicitly document and validate the restriction and direct users to the custom evaluator interface.

**Regression check:** Cover a scalar, a whole-cell struct, and a vector of quadrature-point states for single- and multiple-field values.

## BUG-015 — Threaded-assembly how-to produces a NaN residual

**Severity:** Medium. **Evidence:** Reproduced by executing the how-to.

**Location:** `docs/src/literate_howto/threaded_assembly.jl:29–39`; `src/ExampleElements/HeatEquation.jl:31–36`.

**Problem:** The example calls `work!(assembler, buffer)` without `a`, but `StationaryFourier` calculates its residual from the current temperature gradient. Missing DOF values are intentionally NaN. The stiffness is finite, but the residual is NaN, and adding a body load afterward cannot repair it.

**Failure scenario:** Execute the file. `all(isfinite, K.nzval)` is true; `all(isfinite, r)` is false. Reusing this example to solve a problem yields an unusable right-hand side.

**Smallest reasonable correction:** Initialize `a = zeros(ndofs(dh))` and pass `a` to `work!`.

**Regression check:** Assert finite matrix and residual entries and compare the resulting residual with the explicitly applied body load at the zero initial solution.

## BUG-016 — Empty cell domains fail with internal indexing errors

**Severity:** Low. **Evidence:** Reproduced.

**Location:** `src/states.jl:150–155` (`create_states`); `src/setup.jl:117–119` (`setup_itembuffer`); `src/setup.jl:61` for an empty domain dictionary.

**Problem:** `DomainSpec` intersects the supplied set with its subdomain, which can produce an empty set. `create_states` immediately takes `first(cellset)`; setup later assumes the state dictionary has a first value as well. An entirely empty `setup_domainbuffers` also fails when validation takes its first domain. These failures expose implementation details rather than a clear supported-input contract.

**Reproduction:** `setup_domainbuffer(DomainSpec(dh, nothing, cv; set=Int[]))` throws `BoundsError`. The same can arise from a nonempty set whose intersection with a subdomain is empty.

**Smallest reasonable correction:** Either support empty domains as no-ops, obtaining a buffer prototype without assuming an item exists, or explicitly reject them with an informative `ArgumentError` during setup. Empty multidomain collections should be handled consistently with BUG-008.

**Regression check:** Test explicitly empty sets, empty intersections, and empty domain dictionaries; verify the documented behavior before work starts.

## BUG-017 — Coupling setup incompletely validates domain/task compatibility

**Severity:** Low. **Evidence:** Reproduced for partial domain matches and differing task counts.

**Location:** `src/DomainBuffers.jl:163–172,251–256`; `src/Simulation.jl:66–74`.

**Problem:** The coupling documentation requires matching domains, but setup silently leaves a primary domain completely uncoupled if any partner lacks its key. Work-time routing independently retains whichever partners do have the key, resulting in inconsistent buffer/simulation collections and a delayed `MethodError`. Threaded coupling also indexes partner task-local arrays by primary task number without checking compatibility, yielding a `BoundsError` when the partner has fewer tasks. Neither failure explains the necessary setup correction.

**Failure scenarios:** A and B have `left/right` domains while C has only `left`: coupling A to B and C drops even B's valid `right` coupling, but work still supplies B there. Separately, coupling a four-task primary domain to a one-task partner fails at local index 2.

**Smallest reasonable correction:** Validate matching domain keys and compatible task-local storage at coupling setup and report a descriptive error. If partial overlap is intended to be supported, filter partners consistently in both setup and work routing instead of the current all-or-nothing test. Do not share one mutable partner buffer among concurrent tasks as a shortcut.

**Regression check:** Test partial domain overlap and unequal task counts in both directions; assert either supported, race-free operation or an early informative exception.

## BUG-018 — Public state accessors promise dictionaries but return a limited wrapper

**Severity:** Low. **Evidence:** Reproduced `keys` failure and static documentation mismatch.

**Location:** `src/DomainBuffers.jl:36–61,223–224`; `src/states.jl:1–7`; `docs/src/DomainBuffers/StateVariables.md:36–45`.

**Problem:** `get_state` and `get_old_state` are documented as returning `Dict{Int,S}`, but return `StateVector`. That wrapper implements indexing, assignment, and equality only; it does not support ordinary dictionary operations such as `keys`, `values`, or iteration. The tests themselves reach into `.vals` to traverse states.

**Failure scenario:** Following the documented dictionary contract, `keys(FA.get_state(db))` throws `MethodError`. Postprocessing/checkpoint code written against that contract must use an internal field.

**Smallest reasonable correction:** Document the actual stable state-container API and provide the intended traversal operations. Returning the raw dictionary alone would lose the wrapper's reference-stability benefit when using state flips, so account for that behavior explicitly.

**Regression check:** Exercise the documented iteration/access pattern before and after `update_states!(...; mode=:flip)` and `revert_states!`.

## BUG-019 — Assembler documentation promises nonexistent facet autodiff

**Severity:** Low. **Evidence:** Static dispatch evidence.

**Location:** `docs/src/Workers/Assemblers.md:17–27`; `src/FerriteAssembly.jl:77,107–109`; `src/Autodiff/autodiff.jl`; `src/ItemBuffers/FacetBuffer.jl:57–59`.

**Problem:** The assembler reference says both cell and facet routines fall back to automatic differentiation of their residual functions. Only the cell fallback exists. `facet_routine!` has no default implementation, and requesting an autodiff facet buffer explicitly errors.

**Failure scenario:** A user defines only `facet_residual!` as instructed, then assembles with `KeReAssembler` or a Ferrite assembler. Dispatch fails rather than generating the tangent.

**Smallest reasonable correction:** Correct the reference to state that facet tangents must currently be implemented explicitly. Implementing facet autodiff is a separate feature.

**Regression check:** Keep an explicit test/documentation example for the supported facet interface and verify that unsupported residual-only tangent assembly has a clear diagnostic.

## BUG-020 — WeakForm examples use obsolete loading APIs and the wrong field/sign

**Severity:** Low. **Evidence:** Static API and weak-form comparison.

**Location:** `src/ExampleElements/WeakForm.jl:22–30,42–49`; obsolete package references also appear in `HeatEquation.jl` and `LinearElasticity.jl` docstrings.

**Problem:** Both examples construct `NeumannHandler`, which is not defined/exported by this package; the supported API is `LoadHandler`. They apply body loads to `:c` even though their illustrated unknown is `:u`. The thermal example defines `qn` as outward heat flux but adds a positive Neumann contribution, whereas the displayed thermal weak form has external contribution `-δu*qn`.

**Failure scenario:** Copying the snippets fails at `NeumannHandler`. Renaming the handler alone leaves the body load targeting a nonexistent field (and therefore omitted with a warning). Under the usual convention `r -= fext`, retaining the positive thermal `qn` reverses the intended boundary heat flux.

**Smallest reasonable correction:** Update the snippets to `LoadHandler`, use the actual field name, and show the negative outward thermal-flux contribution with an explicit residual/external-load convention. Replace obsolete `FerriteNeumann.jl` references where the integrated load handler is intended.

**Regression check:** Execute complete versions of the examples and verify both total applied body load and the direction of the thermal boundary contribution.

## BUG-021 — Local-constraint example compares a solution against itself

**Severity:** Low. **Evidence:** Static data-flow proof; the example's existing assertion passes when executed.

**Location:** `docs/src/literate_howto/local_constraints.jl:43–52`.

**Problem:** The purported reference calculation allocates/assembles `K2` and `r2`, then updates `a2` using `K\r`, the original locally constrained system. The final comparison therefore repeats the same solve instead of checking the independent reference. It also calls `apply!` for a residual-increment formulation where the intended comparison should use homogeneous increment constraints.

**Failure scenario:** Break or change the reference `K2/r2` calculation: the final `a2 ≈ a` assertion can still pass because neither reference array participates in the solve.

**Smallest reasonable correction:** Use the independently assembled arrays in the reference solve (`a2 .-= K2 \ r2`) and apply `apply_zero!` to the residual-increment system, matching the first solve's formulation. The introductory claim that `ReAssembler` supports local constraints should also be restricted to `KeReAssembler`.

**Regression check:** Compare the independently computed final solutions and confirm the check fails if the reference stiffness or residual is deliberately perturbed.

## BUG-022 — Some tests labeled as threaded actually run sequentially

**Severity:** Low. **Evidence:** Static setup/dispatch proof.

**Location:** `test/quadpoint_evaluation.jl:26–30`; `test/states.jl:375–377`; `.github/workflows/FerriteMasterCI.yml`.

**Problem:** The single-field quadrature evaluator test loops over `threading` but does not pass it to `setup_domainbuffer`. The state accumulation test says it uses threading but only supplies `colors`; that does not turn threading on. These checks exercise sequential buffers despite their apparent intent. The Ferrite-master workflow also does not explicitly configure multiple Julia threads, unlike the main test job.

**Failure scenario:** A regression affecting single-field threaded quadrature evaluation or the intended threaded state-copy/flip path can escape these particular tests. Other threaded tests exist; this finding does not imply the entire suite lacks parallel coverage.

**Smallest reasonable correction:** Pass the threading flag, assert the resulting buffer type, and exercise state update/revert paths with threaded buffers. Set an explicit multithreaded configuration for the Ferrite-master workflow if that job is intended to cover parallel compatibility.

**Regression check:** Run the corrected tests with at least two worker threads and verify both sequential and threaded cases execute.

## BUG-023 — Figure-factory manifest is incompatible with current package source

**Severity:** Low. **Evidence:** Static project/manifest comparison; no fresh installation attempted.

**Location:** `docs/figure_factory/Manifest.toml:377–395`; `docs/figure_factory/Project.toml`; root `Project.toml`.

**Problem:** The figure-factory manifest resolves Ferrite 0.3.14 and records FerriteAssembly 0.3.3 while pointing to the current repository via `path = "../.."`. Current source is version 0.3.7, requires Ferrite 1, and imports `ConstructionBase`, which is absent from this manifest entry's dependency list. The checked-in environment therefore does not describe a valid dependency graph for the local source it loads.

**Failure scenario:** A developer tries to instantiate/use the figure-factory environment as committed. It supplies an incompatible Ferrite API and incomplete package dependency metadata, requiring manual environment repair before the helper can reliably run.

**Smallest reasonable correction:** Update its project/source configuration and regenerate a compatible manifest, or retire the obsolete environment if the helper is no longer supported. Ensure its visualization dependency can also resolve with Ferrite 1.

**Regression check:** In an isolated clean environment, instantiate the figure-factory project, load the current package, and generate the mixed grid used by the helper.

## Suggested execution order

1. Fix BUG-001–003 first: they affect assembly correctness and can silently corrupt simulations.
2. Fix BUG-004–005 together before regenerating mixed-material tutorial results.
3. Address BUG-006–015 as individual API/runtime repairs, each with its listed regression check.
4. Resolve the low-severity input-validation, documentation, environment, and ineffective-test issues. In particular, do not treat a passing documentation build or the current test suite as verification of the tutorial physics or the uncovered concurrency cases.
