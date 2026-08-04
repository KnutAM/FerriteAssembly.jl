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
    # for the whole cell state (`set_new_to_old_states!` has no default to fall back on).
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
    # `set_new_to_old_states!` calls `copy_state` once per array element (via `map!`),
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
    # overload `copy_state`, to verify that `set_new_to_old_states!` has no fallback and
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
    # `set_new_to_old_states!` must treat it as a single whole-cell-state value (going
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
    # overload, to verify `set_new_to_old_states!` throws `MethodError` for an immutable
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

end

@testset "state variables" begin
    # Defs
    import .TestStateModule: MatA, MatB, MatC, MatD, MatE, MatF, MatG, MatH,
        StateA, StateB, StateC, StateD, StateE, StateF, StateG

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
                for cellnr in 1:getncells(grid)
                    @test states[cellnr] == [StateA(cellnr, i) for i in 1:getnquadpoints(cv)]  # Updated
                end
                states_dc = deepcopy(states) # Allowed to update states during update_states!
                update_states!(container)
                @test old_states == states_dc          # Correctly updated values
                states[1][1] = StateA(0,0)
                @test old_states[1][1] == StateA(1,1)   # But not aliased
                allocs = @allocated update_states!(container)
                @test allocs == 0 # Vector{T} where isbitstype(T) should not allocate (MatA fulfills this)
            end

            # set_new_to_old_states!: states (new) should revert to old_states, old_states unaffected
            for container in (buffer, Simulation(buffer))
                work!(r_assembler, container) # Ensure states holds a known, freshly assembled value
                old_dc = deepcopy(old_states)
                @test states != old_dc # Sanity check that states differ from old_states after work!
                set_new_to_old_states!(container)
                @test states == old_dc          # states reverted to old values
                @test old_states == old_dc      # old_states unaffected
                states[1][1] = StateA(0, 0)
                @test old_states[1][1] == old_dc[1][1] # But not aliased
                allocs = @allocated set_new_to_old_states!(container)
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
            states_dc = deepcopy(states) # Allowed to update states during update_states!
            update_states!(buffer)
            @test old_states == states_dc                  # Correctly updated values
            cellnr = rand(1:getncells(grid))
            coords = getcoordinates(grid, cellnr)
            x_values = [spatial_coordinate(cv, i, coords) for i in 1:getnquadpoints(cv)]
            states[cellnr] = StateB(0, -x_values)
            @test old_states[cellnr] == StateB(cellnr, x_values)   # But not aliased
            allocs = @allocated update_states!(buffer)
            @test allocs == 0 # Vector{T} where !isbitstype(T) should no longer allocate

            # set_new_to_old_states!: states (new) should revert to old_states, old_states unaffected
            # MatB's state is a single mutable struct per cell (not an AbstractArray), so this uses
            # the `copy_state(::StateB) = deepcopy(s)` overload defined above, and allocates.
            old_dc = deepcopy(old_states)
            work!(kr_assembler, buffer)
            @test states != old_dc # Sanity check that states were actually changed by work!
            set_new_to_old_states!(buffer)
            @test states == old_dc          # states reverted to old values
            @test old_states == old_dc      # old_states unaffected
            cellnr = rand(1:getncells(grid))
            states[cellnr].cellnr = -999
            @test old_states[cellnr].cellnr != -999 # But not aliased
            allocs = @allocated set_new_to_old_states!(buffer)
            @test allocs > 0 # Uses the deepcopy-based copy_state overload for non-AbstractArray states

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
            allocs = @allocated update_states!(buffer)
            @test allocs == 0 # Vector{T} where isbitstype(T) should not allocate (MatC fulfills this)

            # set_new_to_old_states!: states (new) should revert to old_states, old_states unaffected
            old_dc = deepcopy(old_states)
            work!(kr_assembler, buffer)
            @test states != old_dc # Sanity check that states were actually changed by work!
            set_new_to_old_states!(buffer)
            @test states == old_dc          # states reverted to old values
            @test old_states == old_dc      # old_states unaffected
            states[1][1] = StateC(999)
            @test old_states[1][1] == old_dc[1][1] # But not aliased
            allocs = @allocated set_new_to_old_states!(buffer)
            @test allocs == 0 # Vector{T} where isbitstype(T) should not allocate (MatC fulfills this)

            # MatE: Vector{StateE} with non-bits elements, exercising the per-element
            # `copy_state` dispatch inside `set_new_to_old_states!`'s `map!` call.
            buffer = setup_domainbuffer(DomainSpec(dh, MatE(), cv))
            states = FerriteAssembly.get_state(buffer)
            old_states = FerriteAssembly.get_old_state(buffer)
            @test isa(old_states, FerriteAssembly.StateVector{Vector{StateE}})
            old_dc = deepcopy(old_states)
            work!(r_assembler, buffer)
            @test states != old_dc # Sanity check that states were actually changed by work!

            TestStateModule.COPY_STATE_ELEM_CALLS[] = 0
            set_new_to_old_states!(buffer)
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
    allocs = @allocated update_states!(buffer)
    @test allocs == 0

    # Smoke-test of set_new_to_old_states! for nothing states (and check no allocations)
    set_new_to_old_states!(buffer) # Compile
    allocs = @allocated set_new_to_old_states!(buffer)
    @test allocs == 0

    gda = DomainSpec(dh, nothing, cv; set=1:getncells(dh.grid)÷2)
    gdb = DomainSpec(dh, nothing, cv; set=setdiff!(Set(1:getncells(dh.grid)), gda.set))
    buffers = setup_domainbuffers(Dict("a"=>gda, "b"=>gdb))
    update_states!(buffers) # Compile
    allocs = @allocated update_states!(buffers)
    @test allocs == 0

    set_new_to_old_states!(buffers) # Compile
    allocs = @allocated set_new_to_old_states!(buffers)
    @test allocs == 0

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
    set_new_to_old_states!(buffer_d)
    @test TestStateModule.COPY_STATE_CALLS[] == getncells(grid_d) # Custom copy_state dispatched for every cell
    @test states_d == old_dc_d          # states reverted to old values
    @test old_states_d == old_dc_d      # old_states unaffected
    cellnr_d = rand(1:getncells(grid_d))
    states_d[cellnr_d].data[1] = -999.0
    @test old_states_d[cellnr_d].data[1] != -999.0 # But not aliased

    # MatF: single mutable struct per cell, like MatB/MatD, but with no `copy_state`
    # overload at all. `set_new_to_old_states!` has no default (deepcopy) fallback to
    # rely on, so it must throw a `MethodError` instead of silently succeeding.
    buffer_f = setup_domainbuffer(DomainSpec(dh_d, MatF(), cv))
    work!(kr_assembler_d, buffer_f)
    @test_throws MethodError set_new_to_old_states!(buffer_f)

    # MatG: cell state is an *immutable* `AbstractVector{StateG}`. Since `ismutable` is
    # false for it, `set_new_to_old_states!` must copy it as a whole (once per cell, via
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
    set_new_to_old_states!(buffer_g)
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
    @test_throws MethodError set_new_to_old_states!(buffer_h)

    # Regression test: a mutable AbstractArray cell state (MatA's Vector{StateA}) whose
    # "new" array has drifted to a different size than the corresponding "old" array must
    # raise a clear `ArgumentError` instead of `map!` silently copying only the common
    # prefix (or erroring obscurely).
    buffer_mismatch = setup_domainbuffer(DomainSpec(dh_d, MatA(), cv))
    states_mismatch = FerriteAssembly.get_state(buffer_mismatch)
    push!(states_mismatch[1], StateA(-1, 0)) # "new" for cell 1 is now longer than "old"
    @test_throws ArgumentError set_new_to_old_states!(buffer_mismatch)
end
