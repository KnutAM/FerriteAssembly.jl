# Could `work_domain_threaded!` use OhMyThreads.jl instead?

This is a design note, not a proposal to add a dependency: CLAUDE.md forbids adding new
`Project.toml` dependencies to this package, and that rule is not being questioned here. The
point of this note is only to record what was actually checked, so the question doesn't need
re-investigating from scratch if the constraint ever changes.

[OhMyThreads.jl](https://github.com/JuliaFolds2/OhMyThreads.jl) (v0.8.9, the latest release at
the time of writing, installed and tested in a scratch environment - **not** added to this
repository's `Project.toml`/`Manifest.toml`) is a general-purpose data-parallelism package:
`tmap`/`tforeach`/`tmapreduce` plus a `@tasks for ... end` macro, over a choice of schedulers
(`DynamicScheduler`, `StaticScheduler`, `GreedyScheduler`, `SerialScheduler`).

Every claim below was verified directly against the installed package (either by reading
`OhMyThreads/src/*.jl` on GitHub, or by running small probe scripts), not reconstructed from
memory - the exact probe scripts are inlined where it matters, so the claims can be re-checked
against a different version.

## Finding 1: `tforeach`/`tmapreduce` with `:dynamic` or `:static` have the *exact same bug*
that `work_domain_threaded!` just had

```julia
using OhMyThreads

function probe()
    ACTIVE = Threads.Atomic{Int}(0)
    SLOWDONE = Threads.Atomic{Bool}(false)
    t0 = time()
    local caught
    try
        tforeach(1:4; scheduler=:static, ntasks=4) do i
            if i == 1
                sleep(0.1)
                error("boom from item 1")
            elseif i == 4
                Threads.atomic_add!(ACTIVE, 1)
                try
                    sleep(1.0); SLOWDONE[] = true
                finally
                    Threads.atomic_sub!(ACTIVE, 1)
                end
            end
        end
    catch e
        caught = e
    end
    return (caught=typeof(caught), elapsed=time()-t0, active=ACTIVE[], slowdone=SLOWDONE[])
end
probe()
# => (caught = TaskFailedException, elapsed ≈ 0.35, active = 1, slowdone = false)
```

`elapsed ≈ 0.35 s` (item 1's delay, not item 4's 1.0 s sleep), `active == 1`: the call returned
while item 4's task was still running. Reading `OhMyThreads/src/implementation.jl` explains
why: both `DynamicScheduler` and `StaticScheduler` spawn one task per chunk and then do
`mapreduce(fetch, op, tasks)` directly, with no step that waits for every task first. `fetch`
on a `Task` that failed raises immediately, so if an earlier-iterated task in that `mapreduce`
happens to fail, later tasks in the list are never waited for - the same orphaned-task hazard
`work_domain_threaded!` had with `Base.Experimental.@sync`, for the same underlying reason.
This is not a version-specific quirk; it's inherent to doing `mapreduce(fetch, op, tasks)`
without a pre-drain.

## Finding 2: `GreedyScheduler` does not have this bug

```julia
tforeach(1:4; scheduler=:greedy, ntasks=4) do i
    # ... same body as above ...
end
# => (caught = TaskFailedException, elapsed ≈ 1.21, active = 0, slowdone = true)
```

`elapsed` now matches the full 1.0 s sleep, and `active == 0`/`slowdone == true`: every task
was actually finished before the call returned. `GreedyScheduler`'s implementation
(`_tmapreduce` for `GreedyScheduler`, same file) explicitly drains every spawned task with
`wait_nothrow` (`wait(t; throw=false)`, or `Base._wait` pre-1.12) *before* its own
`mapreduce(fetch, op, tasks)` step - i.e. it already does the "wait for everyone, then extract
the first failure" pattern that `sync_with_abort` (`task_synchronization.jl`) implements by
hand. It also already throws a plain `TaskFailedException`, never a `CompositeException`
(unlike plain `@sync`, which always wraps - see the PR discussion this note originates from).

## Finding 3: persistent per-task buffer selection *is* available on the plain functional API

`OhMyThreads.WithTaskIndex(f)` wraps a mapping function so that it is called as
`f(taskindex, args...)`, with `taskindex` stable in `1:ntasks` for the lifetime of that task -
exactly what `get_local(itembuffers, taskid)` needs. This does **not** require the `@tasks`
macro (an earlier draft of this note incorrectly claimed it did, based on only having found the
macro-only `@task_index`/`@local` form - `WithTaskIndex` is the functional-API equivalent, and
was confirmed present and working in the pinned v0.8.9, not just in the "stable" docs which may
track a newer release):

```julia
using OhMyThreads

chunk_vector = [collect((4(i-1)+1):(4i)) for i in 1:40]  # 40 pre-formed chunks of 4 items
tforeach(OhMyThreads.WithTaskIndex((taskindex, chunk) -> begin
    # taskindex ∈ 1:4 here, stable per task, usable to index a pre-allocated buffer
    for item in chunk
        # ...
    end
end), chunk_vector; scheduler=:greedy, ntasks=4, chunking=false)
```

Confirmed working (distinct task indices actually used, each handling a share of the 40
chunks) and confirmed to preserve Finding 2's exception-safety with `WithTaskIndex` added
(`elapsed ≈ 1.08 s`, `active == 0`, `slowdone == true`, same probe body as above).

`chunking=false` is what makes `GreedyScheduler` treat `chunk_vector` itself as the sequence of
"items" to greedily pull one at a time from a shared channel - each item being one of *our*
pre-formed chunks (a `Vector{Int}` of cell numbers), rather than letting OhMyThreads impose its
own automatic (evenly-sized) chunking. This maps directly onto the `TaskChunks`/`get_chunk`
pattern already in this package (`src/Multithreading/TaskChunks.jl`): a shared, lock-protected
queue that tasks pull whole chunks from one at a time.

## Finding 4: OhMyThreads' own task-local storage does not help here

Its `@local`/`TaskLocalValue` machinery (see the "Thread-Safe Storage" page in its docs)
explicitly does **not** persist across separate `tforeach`/`@tasks` calls - a fresh set of
buffers is allocated every call, by design. `ThreadedDomainBuffer.itembuffer::TaskLocals` is
deliberately the opposite: allocated once at `setup_domainbuffer` time and reused for the
lifetime of the buffer, across an unbounded number of future `work!` calls, specifically to
avoid repeated allocation. Adopting OhMyThreads would not remove the need for `TaskLocals` (or
`ThreadedDomainBuffer` holding onto it); it would only be asked, via `WithTaskIndex`, which of
*our own* pre-existing slots to use - exactly the role `taskid` plays in the current code.

## Finding 5: no scheduler has built-in early-abort

Nothing in `GreedyScheduler` (or any other scheduler) checks "has some other task already
failed" before pulling and processing the next item. `sync_with_abort`'s `abort` flag - checked
by the loop body between chunks and items so siblings stop wasting work once one has failed -
has no equivalent here; it would have to be added back exactly as it is today; a shared
`Threads.Atomic{Bool}` passed into the closure and checked inside it.

## What adopting it would actually look like

```julia
using OhMyThreads: tforeach, GreedyScheduler, WithTaskIndex

function work_domain_threaded!(workers, sim::AbstractSingleDomainThreadedSim)
    scatter!(sim)
    itembuffers = get_itembuffer(sim)
    scatter!(workers)
    num_tasks = get_num_tasks(sim)
    for chunk_vector in get_chunks(sim)
        abort = Threads.Atomic{Bool}(false) # Finding 5: still needed, nothing provides this
        tforeach(WithTaskIndex((taskid, chunk) -> begin
            itembuffer = get_local(itembuffers, taskid)
            worker = get_local(workers, taskid)
            for itemnr in chunk
                abort[] && break
                reinit_buffer!(itembuffer, sim, itemnr)
                work_single!(worker, itembuffer)
            end
        end), chunk_vector; scheduler = GreedyScheduler(ntasks = num_tasks), chunking = false)
        # Note: not wrapped to set `abort[]` on the first failure - GreedyScheduler's own
        # exception handling (Finding 2) already waits for everyone and surfaces the first
        # failure as a plain TaskFailedException, so that part doesn't need `sync_with_abort`
        # at all. `abort` only helps *other, still-running* tasks stop early - which requires
        # a caller-visible failure signal *before* `tforeach` itself returns, i.e. still an
        # explicit flag set from inside the failing closure, e.g. via a `try/catch` around the
        # body above that sets `abort[] = true` before `rethrow()`-ing. Doing that only
        # matters for *this specific `chunk_vector`'s remaining work*, since a later
        # `chunk_vector` (color) in the outer `for` loop is never reached anyway once
        # `tforeach` has thrown.
    end
    gather!(itembuffers)
    gather!(workers)
end
```
(Illustrative - written against the verified v0.8.9 behavior above, not tested against this
package's actual `Simulation`/`TaskLocals`/`get_chunks` types.)

This sketch is clean and would let us delete the *runtime queue* half of `TaskChunks.jl` (the
`TaskChunks` struct and `get_chunk`, replaced by `GreedyScheduler`'s own channel) and the
drain/rethrow half of `sync_with_abort`. It would **not** let us delete:
- `TaskLocals`/`ThreadedDomainBuffer.itembuffer` (Finding 4),
- the chunk *construction* utilities in `TaskChunks.jl` - `split_in_chunks`, `create_chunks`,
  `convert_chunk` - which `setup_domainbuffer` still needs regardless of what consumes their
  output (an earlier draft of this note overstated this as "delete `TaskChunks.jl`"; only the
  queue half of that file is actually redundant),
- the abort-flag concept (Finding 5), just moved from `sync_with_abort` into an explicit
  `try/catch` around the closure passed to `tforeach`.

## Recommendation

Even setting the `Project.toml` rule aside: this is a real, fairly clean simplification for the
"spawn + wait + extract first failure" half of the problem (Findings 2 and 3 in particular are
more capable than this note originally gave them credit for), but it is partial - the abort
flag and the persistent-buffer management stay regardless, and the win is "replace ~15 lines of
already-correct, already-tested code with a new external dependency" rather than "remove a
whole category of custom logic". Combined with the actual constraint (no new dependencies),
the hand-rolled `sync_with_abort` in `task_synchronization.jl` is the better choice for this
package today. Worth revisiting if that constraint changes and `GreedyScheduler`'s behavior is
still as described above (it is not the default scheduler, so it may see less maintenance
attention than `Dynamic`/`StaticScheduler`).

### Other packages briefly considered, not investigated to this depth

- **Polyester.jl**: very low per-iteration overhead, but no built-in exception-safety or
  task-local-storage story suited to this use case (it targets tight numeric loops, not
  long-lived cross-call state) - would likely need the same `sync_with_abort`/`TaskLocals`
  machinery on top regardless, without Finding 2/3's help.
- **FLoops.jl** / **Dagger.jl**: heavier-weight general task-scheduling frameworks; not an
  obvious fit for "a fixed, small number of persistent workers pulling pre-formed chunks from a
  queue", and not investigated further here.
