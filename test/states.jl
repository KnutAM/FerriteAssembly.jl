module TestStateModule
    using Ferrite, FerriteAssembly
    struct MatA end
    struct StateA # bitstype
        cellnr::Int
        quadnr::Int
    end
    FerriteAssembly.create_cell_state(::MatA, cv, args...) = [StateA(-1, 0) for _ in 1:getnquadpoints(cv)]
    function FerriteAssembly.element_residual!(re, states::Vector{StateA}, ae, ::MatA, cv, buffer)
        cellnr = cellid(buffer)
        for i in 1:getnquadpoints(cv)
            states[i] = StateA(cellnr, i)
        end
    end

    struct MatB{dim} end
    mutable struct StateB{dim} # not bitstype
        cellnr::Int
        const quad_coordinates::Vector{Vec{dim,Float64}}
    end
    FerriteAssembly.create_cell_state(::MatB{dim}, cv, args...) where dim = StateB(-1, [zero(Vec{dim}) for i in 1:getnquadpoints(cv)])
    function FerriteAssembly.element_residual!(re, states::StateB, ae, ::MatB, cv, buffer)
        coords = getcoordinates(buffer)
        states.cellnr = cellid(buffer)
        for i in 1:getnquadpoints(cv)
            x = spatial_coordinate(cv, i, coords)
            states.quad_coordinates[i] = x
        end
    end

    struct MatC end
    struct StateC
        counter::Int
    end
    FerriteAssembly.create_cell_state(::MatC, cv, args...) = [StateC(0) for _ in 1:getnquadpoints(cv)]
    function FerriteAssembly.element_residual!(re, states::Vector{StateC}, ae, ::MatC, cv, buffer)
        old_states = FerriteAssembly.get_old_state(buffer)
        for i in 1:getnquadpoints(cv)
            states[i] = StateC(old_states[i].counter + 1)
        end
    end

    # With contained vectors, comparison gives false in all cases even with equal values...
    Base.:(==)(a::StateB, b::StateB) = (a.cellnr==b.cellnr && (mapreduce((ax, bx)->ax==bx, *, a.quad_coordinates, b.quad_coordinates)))

    # StateB is not `isbits` and not an `AbstractArray`, so `copy_state` must be overloaded
    # for the whole cell state (`revert_states!` has no default to fall back on).
    FerriteAssembly.copy_state(s::StateB) = deepcopy(s)

    # MatD: single mutable struct per cell (like MatB), but with a `copy_state` overload
    # that reuses the state's own copy logic instead of plain `deepcopy`.
    struct MatD end
    mutable struct StateD
        cellnr::Int
        const data::Vector{Float64}
    end
    FerriteAssembly.create_cell_state(::MatD, cv, args...) = StateD(-1, zeros(getnquadpoints(cv)))
    function FerriteAssembly.element_residual!(re, state::StateD, ae, ::MatD, cv, buffer)
        state.cellnr = cellid(buffer)
        fill!(state.data, cellid(buffer))
    end
    Base.:(==)(a::StateD, b::StateD) = (a.cellnr == b.cellnr && a.data == b.data)

    const COPY_STATE_CALLS = Ref(0)
    function FerriteAssembly.copy_state(s::StateD)
        COPY_STATE_CALLS[] += 1
        return StateD(s.cellnr, copy(s.data))
    end

    # MatE: cell state is a `Vector{StateE}` (an `AbstractArray`) of non-bits elements, so
    # `revert_states!` calls `copy_state` once per array element (via `map!`),
    # rather than once per cell as for MatD.
    struct MatE end
    mutable struct StateE
        cellnr::Int
        quadnr::Int
        const marker::Vector{Int}
    end
    FerriteAssembly.create_cell_state(::MatE, cv, args...) = [StateE(-1, i, [0]) for i in 1:getnquadpoints(cv)]
    function FerriteAssembly.element_residual!(re, states::Vector{StateE}, ae, ::MatE, cv, buffer)
        cellnr = cellid(buffer)
        for i in 1:getnquadpoints(cv)
            states[i] = StateE(cellnr, i, [cellnr])
        end
    end
    Base.:(==)(a::StateE, b::StateE) = (a.cellnr == b.cellnr && a.quadnr == b.quadnr && a.marker == b.marker)

    const COPY_STATE_ELEM_CALLS = Ref(0)
    function FerriteAssembly.copy_state(s::StateE)
        COPY_STATE_ELEM_CALLS[] += 1
        return StateE(s.cellnr, s.quadnr, copy(s.marker))
    end

    # MatF: single mutable struct per cell, like MatB/MatD, but deliberately does NOT
    # overload `copy_state`, to verify that `revert_states!` has no fallback and
    # throws a `MethodError` instead of silently `deepcopy`-ing the state.
    struct MatF end
    mutable struct StateF
        cellnr::Int
    end
    FerriteAssembly.create_cell_state(::MatF, cv, args...) = StateF(-1)
    function FerriteAssembly.element_residual!(re, state::StateF, ae, ::MatF, cv, buffer)
        state.cellnr = cellid(buffer)
    end

    # StateG: mutable per-quadpoint state, referenced from (not owned by) an *immutable*
    # `AbstractVector` wrapper below. Mutating a `StateG` fetched via `getindex` still
    # mutates the shared object, even though the wrapper itself can't be `setindex!`-ed.
    mutable struct StateG
        cellnr::Int
        quadnr::Int
    end
    Base.:(==)(a::StateG, b::StateG) = (a.cellnr == b.cellnr && a.quadnr == b.quadnr)

    # ImmutableStates: an immutable `AbstractVector{StateG}` (e.g. backed by an `NTuple`,
    # like `StaticArrays.SVector` would be). `ismutable(::ImmutableStates) == false`, so
    # `revert_states!` must treat it as a single whole-cell-state value (going
    # through `copy_state` for the whole array) rather than element-wise via `map!`.
    struct ImmutableStates{N} <: AbstractVector{StateG}
        data::NTuple{N,StateG}
    end
    Base.size(x::ImmutableStates) = (length(x.data),)
    Base.getindex(x::ImmutableStates, i::Int) = x.data[i]
    Base.IndexStyle(::Type{<:ImmutableStates}) = IndexLinear()

    # MatG: has a `copy_state` overload for the whole `ImmutableStates` array.
    struct MatG end
    FerriteAssembly.create_cell_state(::MatG, cv, args...) = ImmutableStates(ntuple(i -> StateG(-1, i), getnquadpoints(cv)))
    function FerriteAssembly.element_residual!(re, states::ImmutableStates, ae, ::MatG, cv, buffer)
        cellnr = cellid(buffer)
        for s in states # mutate the elements in place; `states` itself is never `setindex!`-ed
            s.cellnr = cellnr
        end
    end

    const COPY_STATE_WHOLE_ARRAY_CALLS = Ref(0)
    function FerriteAssembly.copy_state(s::ImmutableStates)
        COPY_STATE_WHOLE_ARRAY_CALLS[] += 1
        return ImmutableStates(map(x -> StateG(x.cellnr, x.quadnr), s.data))
    end

    # MatH: same immutable-array state shape as MatG, but deliberately has NO `copy_state`
    # overload, to verify `revert_states!` throws `MethodError` for an immutable
    # `AbstractArray` cell state just as it does for a non-array one (MatF).
    struct ImmutableStatesNoOverload{N} <: AbstractVector{StateG}
        data::NTuple{N,StateG}
    end
    Base.size(x::ImmutableStatesNoOverload) = (length(x.data),)
    Base.getindex(x::ImmutableStatesNoOverload, i::Int) = x.data[i]
    Base.IndexStyle(::Type{<:ImmutableStatesNoOverload}) = IndexLinear()

    struct MatH end
    FerriteAssembly.create_cell_state(::MatH, cv, args...) = ImmutableStatesNoOverload(ntuple(i -> StateG(-1, i), getnquadpoints(cv)))
    function FerriteAssembly.element_residual!(re, states::ImmutableStatesNoOverload, ae, ::MatH, cv, buffer)
        cellnr = cellid(buffer)
        for s in states
            s.cellnr = cellnr
        end
    end

    # StateI: *immutable* wrapper around a mutable `Vector` payload (the motivating case for
    # `copy_state!`: `ismutable(::StateI) == false`, but `vals` can still be `copyto!`-ed into
    # in place). Defines *both* `copy_state` and `copy_state!` to verify that `copy_state!`
    # takes precedence when both are applicable.
    struct MatI end
    struct StateI
        vals::Vector{Float64}
    end
    FerriteAssembly.create_cell_state(::MatI, cv, args...) = StateI(fill(-1.0, getnquadpoints(cv)))
    function FerriteAssembly.element_residual!(re, state::StateI, ae, ::MatI, cv, buffer)
        fill!(state.vals, Float64(cellid(buffer)))
    end
    Base.:(==)(a::StateI, b::StateI) = a.vals == b.vals

    const COPY_STATE_I_CALLS = Ref(0)
    FerriteAssembly.copy_state(s::StateI) = (COPY_STATE_I_CALLS[] += 1; StateI(copy(s.vals)))

    const COPY_STATE_BANG_I_CALLS = Ref(0)
    function FerriteAssembly.copy_state!(dst::StateI, src::StateI)
        COPY_STATE_BANG_I_CALLS[] += 1
        copyto!(dst.vals, src.vals)
        return nothing
    end

    # MatJ: cell state is a `Vector{StateJ}` (a mutable `AbstractArray`) whose *elements* are
    # themselves immutable wrappers around a mutable `Vector`, exercising `copy_state!` in the
    # per-array-element branch of `_copy_states!` (as opposed to MatI's whole-cell-state
    # branch). Defines only `copy_state!` (no `copy_state`).
    struct MatJ end
    struct StateJ
        vals::Vector{Float64}
    end
    FerriteAssembly.create_cell_state(::MatJ, cv, args...) = [StateJ([-1.0]) for _ in 1:getnquadpoints(cv)]
    function FerriteAssembly.element_residual!(re, states::Vector{StateJ}, ae, ::MatJ, cv, buffer)
        cellnr = cellid(buffer)
        for s in states
            fill!(s.vals, Float64(cellnr))
        end
    end
    Base.:(==)(a::StateJ, b::StateJ) = a.vals == b.vals

    const COPY_STATE_BANG_J_CALLS = Ref(0)
    function FerriteAssembly.copy_state!(dst::StateJ, src::StateJ)
        COPY_STATE_BANG_J_CALLS[] += 1
        copyto!(dst.vals, src.vals)
        return nothing
    end

    # MatK: cell state is `Vector{StateK}(undef, n)` - deliberately left with unassigned
    # elements until assembly writes into them. The "old" state is never assembled into
    # directly, so it is still fully unassigned the first time `update_states!` is called;
    # `_copy_states!` must not read an unassigned destination array element (only
    # `copy_state` is defined here, no `copy_state!`).
    struct MatK end
    mutable struct StateK
        cellnr::Int
        quadnr::Int
    end
    FerriteAssembly.create_cell_state(::MatK, cv, args...) = Vector{StateK}(undef, getnquadpoints(cv))
    function FerriteAssembly.element_residual!(re, states::Vector{StateK}, ae, ::MatK, cv, buffer)
        cellnr = cellid(buffer)
        for i in 1:getnquadpoints(cv)
            states[i] = StateK(cellnr, i)
        end
    end
    Base.:(==)(a::StateK, b::StateK) = (a.cellnr == b.cellnr && a.quadnr == b.quadnr)
    FerriteAssembly.copy_state(s::StateK) = StateK(s.cellnr, s.quadnr)

end

@testset "state variables" begin
    # Defs
    import .TestStateModule: MatA, MatB, MatC, MatD, MatE, MatF, MatG, MatH, MatI, MatJ, MatK,
        StateA, StateB, StateC, StateD, StateE, StateF, StateG, StateI, StateJ, StateK

    # `update_states!` accepts a keyword argument (`mode`). On Julia versions before 1.12,
    # a bare `@allocated update_states!(x)` (or with `mode=...`) written directly at
    # top-level/testset scope measures a small constant keyword-dispatch overhead unrelated
    # to `update_states!`'s own allocations - this goes away only when the measured call
    # happens inside a compiled function. Route all allocation-sensitive `update_states!`
    # calls through these tiny helpers instead of a bare `allocs = @allocated ...`; as with
    # the previous `allocs = @allocated ...` pattern, the method must already be compiled
    # (e.g. by a preceding real call) before the result is meaningful - these helpers make a
    # single call, exactly like the pattern they replace, so call counts (and thus e.g.
    # `mode = :flip`'s toggling of old/new) are unaffected.
    _flip!(container) = update_states!(container; mode=:flip)
    _alloc_update!(container) = @allocated update_states!(container)
    _alloc_flip!(container) = @allocated update_states!(container; mode=:flip)

    for (CT, Dim) in ((Line, 1), (QuadraticTriangle, 2), (Hexahedron, 3))
        @testset "$CT" begin
            grid = generate_grid(CT, ntuple(_->3, Dim))
            ip = geometric_interpolation(CT)
            dh = DofHandler(grid); add!(dh, :u, ip^Dim); close!(dh);

            K = allocate_matrix(dh)
            r = zeros(ndofs(dh));
            kr_assembler = start_assemble(K, r)
            r_assembler = FerriteAssembly.ReAssembler(r)
            a = copy(r)
            RefShape = Ferrite.getrefshape(ip)
            cv = CellValues(QuadratureRule{RefShape}(2), ip^Dim, ip);

            # MatA: 
            # - Check correct values before and after update
            # - Check unaliased old and new after update_states!
            buffer = setup_domainbuffer(DomainSpec(dh, MatA(), cv))
            states = FerriteAssembly.get_state(buffer)
            old_states = FerriteAssembly.get_old_state(buffer)
            @test isa(old_states, FerriteAssembly.StateVector{Vector{StateA}})
            @test old_states == states
            @test old_states[1] == [StateA(-1, 0) for _ in 1:getnquadpoints(cv)]
            work!(r_assembler, buffer)
            @test old_states[1] == [StateA(-1, 0) for _ in 1:getnquadpoints(cv)] # Unchanged
            for container in (buffer, Simulation(buffer))
                work!(r_assembler, container) # Ensure states holds a known, freshly assembled value
                for cellnr in 1:getncells(grid)
                    @test states[cellnr] == [StateA(cellnr, i) for i in 1:getnquadpoints(cv)]  # Updated
                end
                states_dc = deepcopy(states) # `mode = :copy` (default): old := new, new left untouched
                update_states!(container)
                @test old_states == states_dc          # Correctly updated values
                @test states == states_dc              # `states` (new) untouched by default `mode = :copy`
                states[1][1] = StateA(0,0)
                @test old_states[1][1] == StateA(1,1)   # But not aliased
                @test _alloc_update!(container) == 0 # Vector{T} where isbitstype(T) should not allocate (MatA fulfills this)
            end
            @test_throws ArgumentError update_states!(buffer; mode=:bogus)

            # mode = :flip: reference swap (no copy_state requirement, always allocation-free),
            # but the "new" container ends up holding the values from *before* this call.
            work!(r_assembler, buffer)
            old_before_flip = deepcopy(old_states)
            states_before_flip = deepcopy(states)
            _flip!(buffer)
            @test old_states == states_before_flip # old_states now holds what was `states`
            @test states == old_before_flip        # states now holds what was `old_states` (stale)
            @test _alloc_flip!(buffer) == 0

            # revert_states!: states (new) should revert to old_states, old_states unaffected
            for container in (buffer, Simulation(buffer))
                work!(r_assembler, container) # Ensure states holds a known, freshly assembled value
                old_dc = deepcopy(old_states)
                @test states != old_dc # Sanity check that states differ from old_states after work!
                revert_states!(container)
                @test states == old_dc          # states reverted to old values
                @test old_states == old_dc      # old_states unaffected
                states[1][1] = StateA(0, 0)
                @test old_states[1][1] == old_dc[1][1] # But not aliased
                allocs = @allocated revert_states!(container)
                @test allocs == 0 # Vector{T} where isbitstype(T) should not allocate (MatA fulfills this)
            end

            # MatB (not bitstype)
            # - Check correct values before and after update
            # - Check unaliased old and new after update_states!
            buffer = setup_domainbuffer(DomainSpec(dh, MatB{Dim}(), cv))
            states = FerriteAssembly.get_state(buffer)
            old_states = FerriteAssembly.get_old_state(buffer)
            @test isa(old_states, FerriteAssembly.StateVector{StateB{Dim}})
            @test old_states == states
            @test old_states[1] == StateB(-1, [zero(Vec{Dim}) for i in 1:getnquadpoints(cv)])
            work!(kr_assembler, buffer)
            @test old_states[1] == StateB(-1, [zero(Vec{Dim}) for i in 1:getnquadpoints(cv)]) # Unchanged
            for cellnr in 1:getncells(grid)
                coords = getcoordinates(grid, cellnr)
                x_values = [spatial_coordinate(cv, i, coords) for i in 1:getnquadpoints(cv)]
                @test states[cellnr] == StateB(cellnr, x_values)                          # Updated
            end
            states_dc = deepcopy(states) # `mode = :copy` (default): old := new, new left untouched
            update_states!(buffer)
            @test old_states == states_dc                  # Correctly updated values
            @test states == states_dc                      # `states` (new) untouched by default `mode = :copy`
            cellnr = rand(1:getncells(grid))
            coords = getcoordinates(grid, cellnr)
            x_values = [spatial_coordinate(cv, i, coords) for i in 1:getnquadpoints(cv)]
            states[cellnr] = StateB(0, -x_values)
            @test old_states[cellnr] == StateB(cellnr, x_values)   # But not aliased
            allocs = @allocated update_states!(buffer)
            @test allocs > 0 # MatB's state (not an AbstractArray) uses the deepcopy-based copy_state overload

            # mode = :flip: no `copy_state` requirement, always allocation-free
            work!(kr_assembler, buffer)
            old_before_flip = deepcopy(old_states)
            states_before_flip = deepcopy(states)
            _flip!(buffer)
            @test old_states == states_before_flip # old_states now holds what was `states`
            @test states == old_before_flip        # states now holds what was `old_states` (stale)
            @test _alloc_flip!(buffer) == 0

            # revert_states!: states (new) should revert to old_states, old_states unaffected.
            # MatB's state is a single mutable struct per cell (not an AbstractArray), so this
            # uses the `copy_state(::StateB) = deepcopy(s)` overload defined above, and allocates.
            old_dc = deepcopy(old_states)
            work!(kr_assembler, buffer)
            @test states != old_dc # Sanity check that states were actually changed by work!
            revert_states!(buffer)
            @test states == old_dc          # states reverted to old values
            @test old_states == old_dc      # old_states unaffected
            cellnr = rand(1:getncells(grid))
            states[cellnr].cellnr = -999
            @test old_states[cellnr].cellnr != -999 # But not aliased
            allocs = @allocated revert_states!(buffer)
            @test allocs > 0 # Uses the deepcopy-based copy_state overload for non-AbstractArray states

            # set_new_to_old_states! is a deprecated alias for revert_states!, kept for
            # backwards compatibility - check it still works (and warns).
            old_dc = deepcopy(old_states)
            work!(kr_assembler, buffer)
            @test states != old_dc # Sanity check that states were actually changed by work!
            @test_deprecated set_new_to_old_states!(buffer)
            @test states == old_dc # still forwards correctly to revert_states!
            @test old_states == old_dc

            # MatC (accumulation), using threading as well
            colors = create_coloring(grid)
            buffer = setup_domainbuffer(DomainSpec(dh, MatC(), cv; colors=colors))
            states = FerriteAssembly.get_state(buffer)
            old_states = FerriteAssembly.get_old_state(buffer)
            @test isa(old_states, FerriteAssembly.StateVector{Vector{StateC}})
            @test old_states == states
            @test old_states[1][1] == StateC(0)
            work!(kr_assembler, buffer)
            @test old_states[1][1] == StateC(0)
            @test states[1][1] == StateC(1) # Added 1
            update_states!(buffer)
            @test old_states[1][1] == StateC(1) # Updated
            states[1][1] = StateC(0)        # Set states to zero to test aliasing and update in next assembly
            @test old_states[1][1] == StateC(1) # But not aliased
            work!(kr_assembler, buffer)
            @test states[1][1] == StateC(2) # Added 1 from old_states[1][1] (not from states[1][1] which was StateC(0))
            for cellnr in 1:getncells(grid)
                @test states[cellnr][2] == StateC(2) # Check that all are updated
            end
            @test _alloc_update!(buffer) == 0 # Vector{T} where isbitstype(T) should not allocate (MatC fulfills this)

            # revert_states!: states (new) should revert to old_states, old_states unaffected
            old_dc = deepcopy(old_states)
            work!(kr_assembler, buffer)
            @test states != old_dc # Sanity check that states were actually changed by work!
            revert_states!(buffer)
            @test states == old_dc          # states reverted to old values
            @test old_states == old_dc      # old_states unaffected
            states[1][1] = StateC(999)
            @test old_states[1][1] == old_dc[1][1] # But not aliased
            allocs = @allocated revert_states!(buffer)
            @test allocs == 0 # Vector{T} where isbitstype(T) should not allocate (MatC fulfills this)

            # MatE: Vector{StateE} with non-bits elements, exercising the per-element
            # `copy_state` dispatch inside `revert_states!`'s `map!` call.
            buffer = setup_domainbuffer(DomainSpec(dh, MatE(), cv))
            states = FerriteAssembly.get_state(buffer)
            old_states = FerriteAssembly.get_old_state(buffer)
            @test isa(old_states, FerriteAssembly.StateVector{Vector{StateE}})
            old_dc = deepcopy(old_states)
            work!(r_assembler, buffer)
            @test states != old_dc # Sanity check that states were actually changed by work!

            TestStateModule.COPY_STATE_ELEM_CALLS[] = 0
            revert_states!(buffer)
            nqp = getnquadpoints(cv)
            @test TestStateModule.COPY_STATE_ELEM_CALLS[] == getncells(grid) * nqp # copy_state dispatched per array element
            @test states == old_dc          # states reverted to old values
            @test old_states == old_dc      # old_states unaffected
            states[1][1].marker[1] = -999
            @test old_states[1][1].marker[1] != -999 # But not aliased (element-wise copy, not shared)
        end
    end

    ip = Lagrange{RefTriangle,1}()
    dh = DofHandler(generate_grid(Triangle, (2,2))); add!(dh, :u, ip); close!(dh);
    
    # Smoke-test of update_states! for nothing states (and check no allocations)
    cv = CellValues(QuadratureRule{RefTriangle}(2), ip)
    buffer = setup_domainbuffer(DomainSpec(dh, nothing, cv))
    @test isa(FerriteAssembly.get_state(buffer), FerriteAssembly.StateVector{Vector{Nothing}})
    update_states!(buffer) # Compile
    @test _alloc_update!(buffer) == 0
    _flip!(buffer) # Compile
    @test _alloc_flip!(buffer) == 0
    @test_throws ArgumentError update_states!(buffer; mode=:bogus)

    # Smoke-test of revert_states! for nothing states (and check no allocations)
    revert_states!(buffer) # Compile
    allocs = @allocated revert_states!(buffer)
    @test allocs == 0

    # Smoke-test of the deprecated set_new_to_old_states! alias for nothing states
    @test_deprecated set_new_to_old_states!(buffer)

    gda = DomainSpec(dh, nothing, cv; set=1:getncells(dh.grid)÷2)
    gdb = DomainSpec(dh, nothing, cv; set=setdiff!(Set(1:getncells(dh.grid)), gda.set))
    buffers = setup_domainbuffers(Dict("a"=>gda, "b"=>gdb))
    update_states!(buffers) # Compile
    @test _alloc_update!(buffers) == 0

    # `mode` must thread through the multi-domain `Dict` layer and the `Simulation` wrapper
    # too (the layers the reported issue's regression went through), and so must revert_states!
    # and its deprecated `set_new_to_old_states!` alias.
    for container in (buffers, Simulation(buffers))
        _flip!(container) # Compile
        @test _alloc_flip!(container) == 0
        revert_states!(container) # Compile
        allocs = @allocated revert_states!(container)
        @test allocs == 0
        @test_deprecated set_new_to_old_states!(container)
    end

    # MatD: single mutable struct per cell, overloading `FerriteAssembly.copy_state`
    # with logic other than plain `deepcopy` (contrast with MatB above).
    grid_d = generate_grid(Triangle, (2, 2))
    dh_d = DofHandler(grid_d); add!(dh_d, :u, ip); close!(dh_d)
    K_d = allocate_matrix(dh_d)
    r_d = zeros(ndofs(dh_d))
    kr_assembler_d = start_assemble(K_d, r_d)
    buffer_d = setup_domainbuffer(DomainSpec(dh_d, MatD(), cv))
    states_d = FerriteAssembly.get_state(buffer_d)
    old_states_d = FerriteAssembly.get_old_state(buffer_d)
    @test isa(old_states_d, FerriteAssembly.StateVector{StateD})

    old_dc_d = deepcopy(old_states_d)
    work!(kr_assembler_d, buffer_d)
    @test states_d != old_dc_d # Sanity check that states were actually changed by work!

    TestStateModule.COPY_STATE_CALLS[] = 0
    revert_states!(buffer_d)
    @test TestStateModule.COPY_STATE_CALLS[] == getncells(grid_d) # Custom copy_state dispatched for every cell
    @test states_d == old_dc_d          # states reverted to old values
    @test old_states_d == old_dc_d      # old_states unaffected
    cellnr_d = rand(1:getncells(grid_d))
    states_d[cellnr_d].data[1] = -999.0
    @test old_states_d[cellnr_d].data[1] != -999.0 # But not aliased

    # MatF: single mutable struct per cell, like MatB/MatD, but with no `copy_state`
    # overload at all. `revert_states!` has no default (deepcopy) fallback to
    # rely on, so it must throw a `MethodError` instead of silently succeeding.
    buffer_f = setup_domainbuffer(DomainSpec(dh_d, MatF(), cv))
    work!(kr_assembler_d, buffer_f)
    @test_throws MethodError revert_states!(buffer_f)
    # ...and therefore also the default `mode = :copy` of `update_states!`; `mode = :flip`
    # has no `copy_state` requirement and keeps working for such a state type.
    @test_throws MethodError update_states!(buffer_f)
    update_states!(buffer_f; mode=:flip)

    # MatG: cell state is an *immutable* `AbstractVector{StateG}`. Since `ismutable` is
    # false for it, `revert_states!` must copy it as a whole (once per cell, via
    # the `copy_state(::ImmutableStates)` overload above) rather than element-wise via
    # `map!` (which requires a mutable destination).
    buffer_g = setup_domainbuffer(DomainSpec(dh_d, MatG(), cv))
    states_g = FerriteAssembly.get_state(buffer_g)
    old_states_g = FerriteAssembly.get_old_state(buffer_g)
    @test !ismutable(old_states_g[1]) # Sanity check: this is the scenario being tested

    old_dc_g = deepcopy(old_states_g)
    work!(kr_assembler_d, buffer_g)
    @test states_g != old_dc_g # Sanity check that states were actually changed by work!

    TestStateModule.COPY_STATE_WHOLE_ARRAY_CALLS[] = 0
    revert_states!(buffer_g)
    @test TestStateModule.COPY_STATE_WHOLE_ARRAY_CALLS[] == getncells(grid_d) # Whole array copied once per cell, not per element
    @test states_g == old_dc_g          # states reverted to old values
    @test old_states_g == old_dc_g      # old_states unaffected
    states_g[1][1].cellnr = -999
    @test old_states_g[1][1].cellnr != -999 # But not aliased

    # MatH: same immutable-array cell state shape as MatG, but with no `copy_state`
    # overload, to verify that an immutable `AbstractArray` cell state without an overload
    # throws `MethodError` just like a non-array one (MatF), rather than erroring inside
    # `map!` from trying to mutate an immutable destination.
    buffer_h = setup_domainbuffer(DomainSpec(dh_d, MatH(), cv))
    work!(kr_assembler_d, buffer_h)
    @test_throws MethodError revert_states!(buffer_h)

    # MatI: whole-cell-state is an *immutable* wrapper around a mutable `Vector` payload
    # (`ismutable(::StateI) == false`, but `copy_state!` can still `copyto!` into `vals` in
    # place). Defines both `copy_state` and `copy_state!` - only the latter should be used.
    buffer_i = setup_domainbuffer(DomainSpec(dh_d, MatI(), cv))
    states_i = FerriteAssembly.get_state(buffer_i)
    old_states_i = FerriteAssembly.get_old_state(buffer_i)
    work!(kr_assembler_d, buffer_i)
    @test states_i != old_states_i # Sanity check that states were actually changed by work!

    old_vals_i = old_states_i[1].vals # capture the inner Vector's identity before copying
    TestStateModule.COPY_STATE_I_CALLS[] = 0
    TestStateModule.COPY_STATE_BANG_I_CALLS[] = 0
    update_states!(buffer_i) # mode = :copy (default): old := new
    @test TestStateModule.COPY_STATE_BANG_I_CALLS[] == getncells(grid_d) # copy_state! dispatched once per cell
    @test TestStateModule.COPY_STATE_I_CALLS[] == 0                     # copy_state not used: copy_state! takes precedence
    @test old_states_i == states_i          # old_states updated to the (just-converged) states
    @test old_states_i[1].vals === old_vals_i # not reallocated: same Vector object, mutated in place
    @test _alloc_update!(buffer_i) == 0     # in-place copyto! (no copy_state allocation) is allocation-free

    # revert_states! (the opposite direction) must also use copy_state! in preference to
    # copy_state, and stay allocation-free.
    new_vals_i = states_i[1].vals
    fill!(new_vals_i, -999.0) # corrupt "new" (in place) so revert_states! has something to fix
    TestStateModule.COPY_STATE_I_CALLS[] = 0
    TestStateModule.COPY_STATE_BANG_I_CALLS[] = 0
    revert_states!(buffer_i)
    @test TestStateModule.COPY_STATE_BANG_I_CALLS[] == getncells(grid_d)
    @test TestStateModule.COPY_STATE_I_CALLS[] == 0
    @test states_i == old_states_i
    @test states_i[1].vals === new_vals_i # not reallocated: same Vector object, mutated in place
    allocs = @allocated revert_states!(buffer_i)
    @test allocs == 0

    # MatJ: cell state is a `Vector{StateJ}` (mutable `AbstractArray`) of elements that are
    # themselves immutable wrappers around a mutable `Vector`, exercising `copy_state!` in the
    # per-array-element branch instead of MatI's whole-cell-state branch. Defines only
    # `copy_state!` (no `copy_state` fallback).
    buffer_j = setup_domainbuffer(DomainSpec(dh_d, MatJ(), cv))
    states_j = FerriteAssembly.get_state(buffer_j)
    old_states_j = FerriteAssembly.get_old_state(buffer_j)
    work!(kr_assembler_d, buffer_j)
    @test states_j != old_states_j # Sanity check that states were actually changed by work!

    old_vals_j = old_states_j[1][1].vals # capture the inner Vector's identity before copying
    nqp_j = getnquadpoints(cv)
    TestStateModule.COPY_STATE_BANG_J_CALLS[] = 0
    update_states!(buffer_j) # mode = :copy (default): old := new
    @test TestStateModule.COPY_STATE_BANG_J_CALLS[] == getncells(grid_d) * nqp_j # copy_state! dispatched per array element
    @test old_states_j == states_j          # old_states updated to the (just-converged) states
    @test old_states_j[1][1].vals === old_vals_j # not reallocated: same Vector object, mutated in place
    @test _alloc_update!(buffer_j) == 0     # in-place copyto! is allocation-free

    # revert_states! (the opposite direction) exercises copy_state! per array element too.
    new_vals_j = states_j[1][1].vals
    fill!(new_vals_j, -999.0) # corrupt "new" (in place) so revert_states! has something to fix
    TestStateModule.COPY_STATE_BANG_J_CALLS[] = 0
    revert_states!(buffer_j)
    @test TestStateModule.COPY_STATE_BANG_J_CALLS[] == getncells(grid_d) * nqp_j
    @test states_j == old_states_j
    @test states_j[1][1].vals === new_vals_j # not reallocated: same Vector object, mutated in place
    allocs = @allocated revert_states!(buffer_j)
    @test allocs == 0

    # MatK: regression test for an initially-unassigned destination array element (the "old"
    # state is never assembled into, so its Vector{StateK}(undef, n) elements are genuinely
    # unassigned before the first update_states! call) - must not throw UndefRefError.
    buffer_k = setup_domainbuffer(DomainSpec(dh_d, MatK(), cv))
    old_states_k = FerriteAssembly.get_old_state(buffer_k)
    @test !isassigned(old_states_k[1], 1) # sanity check: this is the scenario being tested
    work!(kr_assembler_d, buffer_k)
    update_states!(buffer_k) # must not throw UndefRefError reading the unassigned destination
    @test isassigned(old_states_k[1], 1)
    @test old_states_k == FerriteAssembly.get_state(buffer_k)

    # Regression test: a mutable AbstractArray cell state (MatA's Vector{StateA}) whose
    # "new" array has drifted to a different size than the corresponding "old" array must
    # raise a clear `ArgumentError` instead of `map!` silently copying only the common
    # prefix (or erroring obscurely).
    buffer_mismatch = setup_domainbuffer(DomainSpec(dh_d, MatA(), cv))
    states_mismatch = FerriteAssembly.get_state(buffer_mismatch)
    push!(states_mismatch[1], StateA(-1, 0)) # "new" for cell 1 is now longer than "old"
    @test_throws ArgumentError revert_states!(buffer_mismatch)
    @test_throws ArgumentError update_states!(buffer_mismatch) # default mode = :copy
    update_states!(buffer_mismatch; mode=:flip) # mode = :flip never touches individual elements
end
