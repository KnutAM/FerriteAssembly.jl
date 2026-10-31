# PR110 review: `fix/threaded-work-wait-for-all-tasks`

Reviewed with a dual-review process: an independent read-only Codex reviewer checked the
review plan and the resulting diff (via the `$independent-review` skill); every finding was
independently re-verified in this session before being accepted.

## Summary

The branch fixes a real race in `work_domain_threaded!` (`src/work.jl`): the old
`Base.Experimental.@sync` rethrew the first failing spawned task's exception without waiting
for its siblings to finish, so those siblings kept mutating shared, persistent task-local
buffers (and `K`/`r`) concurrently with a caller's retry. It replaces `@sync` with a new
primitive, `spawn_tasks_safe_wait`/`drain_tasks`
(`src/Multithreading/task_synchronization.jl`), which always drains every spawned task before
returning or rethrowing.

This review found and fixed **one more real, confirmed bug** in `drain_tasks` (below). With
that fix in place, the branch is correct as far as this review could establish, and the
fix adds no measurable steady-state overhead relative to `main`.

## Correctness review

### Bug found and fixed: `drain_tasks` could stop waiting on a still-running task

**Location:** `src/Multithreading/task_synchronization.jl`, `drain_tasks` (as it was before
this review's fix).

**The code:**

```julia
for t in tasks
    while true
        try
            wait(t)
        catch e
            if !isa(e, TaskFailedException)
                push!(async_exceptions, e)
                continue
            end
        end
        break
    end
end
```

**The bug:** the loop used `isa(e, TaskFailedException)` as a proxy for "`wait(t)` discovered
that `t` itself failed, so move on to the next task." But a `TaskFailedException` can also
arrive *asynchronously* at the waiting task - exactly like the `InterruptException` case the
code already explicitly guards against (and that is exercised by the existing "InterruptException
delivered to the waiting caller" test) - except carrying a `TaskFailedException` that wraps a
*different*, unrelated task. The old code could not tell the two cases apart: on catching such
an exception it took the `break` path, silently discarding the exception (it is not even
collected into `async_exceptions`) and moving on as if `t` had been dealt with, while `t` could
still be running.

**Concrete reproduction** (confirmed independently in this session, matching Codex's
finding): with a child task gated on a `release` flag, and a second, unrelated task made to
schedule `TaskFailedException(unrelated_task)` into the waiting caller while the child is
still running, `drain_tasks` returned successfully with the child still unfinished
(`finished[] == false` immediately after return). This is the same class of bug the whole
branch exists to fix, just one level deeper - a task left running past the point where the
caller believes everything has been drained.

**How realistic this is for `work_domain_threaded!`'s own call site:** Julia itself never
originates an async `TaskFailedException` out of nothing; something has to `schedule(task,
TaskFailedException(x); error = true)`. The only thing Base does that way is `Ctrl-C`
(`InterruptException`), which this code already handles correctly. A `TaskFailedException`
delivered asynchronously requires some *other*, external supervising/cancellation code
forwarding a failure from one task into another. That isn't exercised anywhere in this
package today, so it is not a bug that bites from ordinary use of `work!`. It is, however, a
real type-confusion gap in a primitive explicitly designed to be robust against asynchronous
exceptions in general (its own docstring already distinguished "did it come from waiting on
`t`" from "did it come from something else" - the implementation just didn't check the one
field that actually proves which case it is), so it is worth closing rather than leaving as a
latent trap for the next caller or the next kind of asynchronous exception.

**Fix applied:** check `e.task === t` in addition to the type, so only a `TaskFailedException`
that actually wraps the task currently being waited on is treated as "this task has failed";
anything else (including a `TaskFailedException` for a different task) is collected into
`async_exceptions` and the wait on `t` is retried, exactly like any other asynchronous
exception.

**Regression test added:** `test/task_synchronization.jl`, testset "an unrelated
TaskFailedException delivered asynchronously does not skip draining" - mirrors the existing
`InterruptException`-to-the-caller test exactly, but delivers a `TaskFailedException` wrapping
an unrelated, already-failed task instead, and asserts both that the slow sibling had actually
finished by the time `spawn_tasks_safe_wait` returned, and that the exception's `.task` is
correctly preserved as the unrelated task (not confused with the real failure).

### Everything else checked, no further issues found

1. **CAS-before-`abort[]` ordering** (`spawn_tasks_safe_wait`'s catch block): the CAS on
   `first_failed` happens strictly before `abort[] = true`, on the same task, so any sibling
   that observes `abort[] == true` is guaranteed the CAS has already completed. This is the
   fix from the second round of review in this branch's own history (`.claude/tasks.md`) and
   still holds.
2. **`drain_tasks` otherwise waits for every task exactly once**, via a retry loop that only
   exits on a genuine success or a genuine failure of that specific task (after the fix
   above); asynchronous exceptions are always collected and waiting resumes.
3. **Failure-vs-async-exception precedence** in `spawn_tasks_safe_wait`: a task's own failure
   (via `first_failed[]`) always takes precedence over anything in `async_exceptions`, and the
   thrown exception is always a plain `TaskFailedException`, never a `CompositeException`.
4. **Abort-check granularity** in `work_domain_threaded!` (`src/work.jl`): `abort[]` is checked
   between items and between chunks, never in the middle of a single item's
   `reinit_buffer!`/`work_single!` pair, so no task is ever interrupted mid-item leaving a
   persistent task-local buffer half-written.
5. **`gather!(itembuffers)`/`gather!(workers)`** are only reached if every `chunk_vector`'s
   `spawn_tasks_safe_wait` call returned normally; any throw skips them, so no gather runs
   against buffers that other tasks might still be touching.
6. **No other call site** has the same pattern: `Base.Experimental.@sync` combined with
   `Threads.@spawn` only ever appeared at this one place in `src/`
   (`work_domain_threaded!`); it is now gone entirely.
7. **No downstream code matches on exception type**: nothing in `src/` or `test/` catches
   `CompositeException` from `work!`/`work_domain_threaded!`, so swapping `@sync`'s
   `CompositeException` for a plain `TaskFailedException` does not break any existing caller.
8. **`TaskChunks`** (`src/Multithreading/TaskChunks.jl`, unchanged by this branch) is correctly
   lock-protected and unaffected by the synchronization change.

### Test execution

- `Pkg.test()`: full suite passes, including the new regression test
  (`spawn_tasks_safe_wait`: 17/17, `Threaded work! finishes all tasks before propagating a
  failure (integration)`: 8/8).
- Stress run: full suite run standalone at `-t 1`, `-t 2`, `-t 4`, `-t 12` (this machine has 12
  cores), plus 3 additional repeats each at `-t 4` and `-t 12` (10 full-suite runs total) - all
  clean, no flakiness observed.
- `docs/make.jl`: full docs build (all tutorials, how-tos, and doctests) succeeds with no
  errors.

## Benchmarks vs `main`

**Goal:** quantify whether replacing `Base.Experimental.@sync` with the hand-rolled
`spawn_tasks_safe_wait` (a `Vector{Task}` plus two atomics and a manual drain loop) adds any
steady-state overhead to the common, non-failing threaded-assembly path.

**Method:** a throwaway script (not added to the package or its `Project.toml`), building a
`79,059`-dof / `23,328`-cell hexahedral linear-elasticity problem, timing `work!` with
`@elapsed` (one warmup call to compile, then 5 timed repeats, median reported), for
`num_tasks` in `{1, 2, 4, 8, 12}`. Run identically on this branch and then on `main` (switched
via `git checkout` after stashing this branch's uncommitted fix; no files were lost - the
stash was popped back afterwards), both at `-t 12` on the same machine.

| num_tasks | branch median [s] | main median [s] | difference |
|-----------|-------------------|------------------|------------|
| 1         | 0.624             | 0.629            | ~1% faster |
| 2         | 0.365             | 0.350            | ~4% slower |
| 4         | 0.204             | 0.206            | ~1% faster |
| 8         | 0.129             | 0.122            | ~6% slower |
| 12        | 0.121             | 0.119            | ~2% slower |

**Conclusion:** the differences are within normal run-to-run scheduling noise for this
machine/workload (repeat runs of the same script on the same branch vary by a similar amount -
see the raw 5-sample timings each median was taken from, printed by the script). There is no
measurable steady-state overhead from `spawn_tasks_safe_wait` relative to
`Base.Experimental.@sync`: both scale similarly with `num_tasks`, and neither consistently beats
the other across the five task counts tested.

## Codex review findings

- **Plan review (before implementation):** 1 high-severity finding - the `drain_tasks`
  async-`TaskFailedException` bug described above. **Accepted and fixed**, with an independent
  reproduction in this session confirming both the bug and the fix.
- **Final diff review (after implementation):** no further findings (`NO_FINDINGS`).

## Remaining risks

- The benchmark compares one problem size/shape (a large hex mesh, linear elasticity) on one
  12-core machine; it does not rule out overhead on machines with much higher core counts or
  much smaller per-item workloads (where the per-chunk `Vector{Task}` allocation and atomic
  operations could matter more relative to the work being abort-checked). Given the measured
  parity here and the small, fixed-size nature of the added work (one `Vector{Task}` of length
  `num_tasks` and a couple of atomics per `work!` call, not per item), this is judged unlikely
  to matter in practice.
- The fixed `drain_tasks` bug's realistic trigger (an asynchronously delivered
  `TaskFailedException` referencing an unrelated task) does not occur anywhere in this
  package's own code today; it was a latent robustness gap rather than an observed failure.
  The fix and its regression test close it regardless.
