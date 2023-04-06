@testset "state variables" begin
    # Defs
    struct MatStateTest{T}
        A::Matrix{T}
    end
    struct StateTest
        c::Int
    end
    FA.create_cell_state(::MatStateTest, cv, args...) = [StateTest(0) for _ in 1:getnquadpoints(cv)]
    function FA.element_residual!(re, states::Vector{StateTest}, ae, m::MatStateTest, cv, args...)
        re .= m.A*ae
        for i in 1:getnquadpoints(cv)
            states[i] = StateTest(states[i].c + 1)
        end
    end
    for (CT, Dim) in ((Line, 1), (QuadraticTriangle, 2), (Hexahedron, 3))
        grid = generate_grid(CT, ntuple(_->1, Dim))
        dh = DofHandler(grid); add!(dh, :u, Dim); close!(dh);
        ndpc = Ferrite.ndofs_per_cell(dh); 
        mat=MatStateTest(rand(ndpc,ndpc));
        K = create_sparsity_pattern(dh)
        r = zeros(ndofs(dh));
        a = copy(r)
        ip = Ferrite.default_interpolation(getcelltype(grid))
        RefShape = Ferrite.getrefshape(ip)
        cv = CellVectorValues(QuadratureRule{Dim,RefShape}(2), ip);
        states = FA.create_states(dh, mat, cv);
        cellbuffer = FA.setup_cellbuffer(dh, cv, mat);
        assembler = start_assemble(K,r)
        FA.doassemble!(assembler, cellbuffer, states, dh, a);
        @show ndpc
        @test first(first(states)).c == 1
    end
end