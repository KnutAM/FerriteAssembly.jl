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

## Finding 2: `GreedyScheduler` does not have *that* bug - but has two narrower ones of its own

```julia
tforeach(1:4; scheduler=:greedy, ntasks=4) do i
    # ... same body as above ...
end
# => (caught = TaskFailedException, elapsed ≈ 1.21, active = 0, slowdone = true)
```

`elapsed` now matches the full 1.0 s sleep, and `active == 0`/`slowdone == true`: every task
was actually finished before the call returned when an *ordinary child task* fails.
`GreedyScheduler`'s implementation (`_tmapreduce` for `GreedyScheduler`, same file) explicitly
drains every spawned task with `wait_nothrow` (`wait(t; throw=false)`, or `Base._wait`
pre-1.12) *before* its own `mapreduce(fetch, op, tasks)` step, and already throws a plain
`TaskFailedException`, never a `CompositeException` (unlike plain `@sync`, which always wraps).

However, `wait_nothrow` only suppresses a task's *own* failure being rethrown by `wait` - it
does nothing about an `InterruptException` delivered to the *waiting caller itself* (e.g.
Ctrl-C), which is a different mechanism entirely. Repeating Finding 1's probe but interrupting
the caller instead of making item 1 fail (and gating item 4's "slow" branch on an explicit
`STARTED` flag so the interrupt reliably fires while it is actually running, not before):

```julia
tforeach(1:4; scheduler=:greedy, ntasks=4) do i
    if i == 4
        STARTED[] = true
        sleep(1.0); SLOWDONE[] = true
    end
end
# (after firing an InterruptException at the waiting caller once STARTED[] is seen)
# => (caught = InterruptException, elapsed ≈ 0.0001, active_immediately = 1, slowdone = false)
```

Same orphaned-task symptom as Finding 1, just triggered by a different cause:
`sync_with_abort`'s interrupt-vs-child-failure handling (added specifically in response to a
review finding on this same PR) is not something `GreedyScheduler` provides for free.

Second, when *multiple* tasks fail, `GreedyScheduler` does not reliably report the one that
actually failed first. Its final `mapreduce(fetch, op, tasks)` iterates the spawned *worker*
tasks in their spawn-index order, not in the order their *items* failed - and under
`chunking=false` greedy pulling, which worker ends up processing which item is not
predictable. Repeating `sync_with_abort`'s own "item 3 fails, item 1 waits for that and fails
0.1 s later" scenario (test/task_synchronization.jl) eight times against
`tforeach(...; scheduler=:greedy, ntasks=4, chunking=false)`:

```
1: boom from 1   2: boom from 3   3: boom from 1   4: boom from 1
5: boom from 1   6: boom from 1   7: boom from 1   8: boom from 1
```

Item 3 fails first by construction every time, yet it was reported in only 1 of 8 trials.
`sync_with_abort`'s CAS-based `first_failed` (added for exactly this reason, after an earlier
review finding on this same PR flagged a naive `findfirst`-style approach as equally
unreliable) does not have this problem.

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
chunks) and confirmed to preserve Finding 2's ordinary-child-failure draining with
`WithTaskIndex` added (`elapsed ≈ 1.08 s`, `active == 0`, `slowdone == true`, same probe body
as above) - Finding 2's two caveats (caller-interrupt safety, failure-order under multiple
failures) are unaffected by `WithTaskIndex` either way; they're properties of the scheduler,
not of how the mapping function receives its task index.

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
        end), chunk_vector; scheduler = GreedyScheduler(ntasks = num_tasks, chunking = false))
        # (`chunking` must be passed to the `GreedyScheduler(...)` constructor, not as a
        # separate `tforeach` keyword - OhMyThreads rejects a plain-instance `scheduler=` mixed
        # with other keywords, confirmed: `ArgumentError: Providing an explicit scheduler as
        # well as direct keyword arguments ... is currently not supported`.)
        #
        # Note: `abort[]` here only helps *other, still-running* tasks stop early once one has
        # failed - it does not, by itself, fix Finding 2's two caveats. Getting a caller
        # InterruptException to actually abort cleanly, and getting the *first actual* failure
        # (not just whichever worker's `fetch` happens to run first) reported when more than
        # one task fails, would still need sync_with_abort's drain-with-retry loop and CAS
        # wrapped around this `tforeach` call - at which point most of what this sketch saves
        # relative to Part A has been added back.
    end
    gather!(itembuffers)
    gather!(workers)
end
```
(Illustrative - written against the verified v0.8.9 behavior above, not tested against this
package's actual `Simulation`/`TaskLocals`/`get_chunks` types.)

This sketch is clean and would let us delete the *runtime queue* half of `TaskChunks.jl` (the
`TaskChunks` struct and `get_chunk`, replaced by `GreedyScheduler`'s own channel) and the
*ordinary-failure* half of `sync_with_abort`'s drain/rethrow logic (an earlier draft of this
note claimed the *whole* drain/rethrow half, which Finding 2's two caveats rule out). It would
**not** let us delete:
- `TaskLocals`/`ThreadedDomainBuffer.itembuffer` (Finding 4),
- the chunk *construction* utilities in `TaskChunks.jl` - `split_in_chunks`, `create_chunks`,
  `convert_chunk` - which `setup_domainbuffer` still needs regardless of what consumes their
  output (an earlier draft of this note overstated this as "delete `TaskChunks.jl`"; only the
  queue half of that file is actually redundant),
- the abort-flag concept (Finding 5), just moved from `sync_with_abort` into an explicit
  `try/catch` around the closure passed to `tforeach`,
- the caller-interrupt retry loop or the CAS-based true-first-failure tracking (Finding 2) -
  `sync_with_abort`'s drain loop would still need to wrap the `tforeach` call to get either of
  those back, exactly as sketched above.

## Recommendation

Even setting the `Project.toml` rule aside: this is a real but narrower simplification than it
first looked. `GreedyScheduler` genuinely removes the "wait for everyone before an *ordinary*
failure propagates" part by itself (Finding 2), and `WithTaskIndex` genuinely removes the need
for the `@tasks` macro for persistent-buffer indexing (Finding 3) - those two corrections to
this note's own earlier drafts are real. But the interrupt-safety and true-first-failure
guarantees `sync_with_abort` provides (both added in direct response to review findings on this
same PR) are not things `GreedyScheduler` provides, so adopting OhMyThreads would not let this
package delete the part of `sync_with_abort` that those guarantees live in - only the "ordinary
wait" part, with the chunking/queue machinery as the main additional win. Combined with the
actual constraint (no new dependencies) and that `GreedyScheduler` is not OhMyThreads' default
scheduler (so may see less maintenance attention than `Dynamic`/`StaticScheduler`, which do
have Finding 1's bug), the hand-rolled `sync_with_abort` in `task_synchronization.jl` remains
the better choice for this package today.

### Other packages briefly considered, not investigated to this depth

- **Polyester.jl**: very low per-iteration overhead, but no built-in exception-safety or
  task-local-storage story suited to this use case (it targets tight numeric loops, not
  long-lived cross-call state) - would likely need the same `sync_with_abort`/`TaskLocals`
  machinery on top regardless, without Finding 2/3's help.
- **FLoops.jl** / **Dagger.jl**: heavier-weight general task-scheduling frameworks; not an
  obvious fit for "a fixed, small number of persistent workers pulling pre-formed chunks from a
  queue", and not investigated further here.
