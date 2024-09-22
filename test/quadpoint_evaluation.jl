@testset "QuadPointEvaluator" begin
    function dofhandler(grid, fields::Pair{Symbol, <:Interpolation}...)
        dh = DofHandler(grid)
        for (key, ip) in fields 
            add!(dh, key, ip)
        end
        return close!(dh)
    end
    # Setup a grid and dofhandler for testing
    # "left" and "right" cellsets to test dual domains
    grid = generate_grid(Quadrilateral, (4,3))
    addcellset!(grid, "left", x -> x[1] < 1e-4)
    addcellset!(grid, "right", setdiff(1:getncells(grid), getcellset(grid, "left")))
    struct QEMat{K} end # Dummy type to dispatch for evaluation
    FerriteAssembly.create_cell_state(m::QEMat, cv::AbstractCellValues, args...) = [qe_state(m) for _ in 1:getnquadpoints(cv)]
    FerriteAssembly.create_cell_state(m::QEMat, cv::NamedTuple, args...) = FerriteAssembly.create_cell_state(m, first(cv))
    qe_state(::QEMat) = nothing

    ipu = Lagrange{RefQuadrilateral, 1}()
    ipv = ipu^2
    dh1 = dofhandler(grid, :u => ipu)
    dh2 = dofhandler(grid, :u => ipu, :v => ipv)
    qr = QuadratureRule{RefQuadrilateral}(2)
    cvu = CellValues(qr, ipu)
    cvv = CellValues(qr, ipv)
    
    @testset "SingleDomain" begin
        @testset "Single field" begin
            qe_state(::QEMat{1}) = rand()
            db = setup_domainbuffer(DomainSpec(dh1, QEMat{1}(), cvu))
            states = FerriteAssembly.get_state(db)
            @assert states[1][1] != states[2][1] # Catch bugs if all cells would be created equal. 
    
            foo(m, u, ∇u, qp_state) = 2 * qp_state
            qe = QuadPointEvaluator{Float64}(db, foo)
            work!(qe, db)
            for i in 1:getncells(grid)
                @test qe.data[i] ≈ 2 * states[i]
            end
            bar(m, u, ∇u, qp_state) = u
            a = rand(ndofs(dh1))
            qe_f = QuadPointEvaluator{Float64}(db, bar)
            work!(qe_f, db; a)
            proj = L2Projector(ipu, grid)
            ap = project(proj, qe_f.data, qr)

            @test evaluate_at_grid_nodes(dh1, a, :u) ≈ evaluate_at_grid_nodes(proj, ap)
        end
        @testset "Multiple fields" begin
            qe_state(::QEMat{2}) = rand()
            db = setup_domainbuffer(DomainSpec(dh2, QEMat{2}(), (u = cvu, v = cvv)))
            states = FerriteAssembly.get_state(db)
    
            foo(m, u, ∇u, qp_state) = 2 * qp_state
            qe = QuadPointEvaluator{Float64}(db, foo)
            work!(qe, db)
            for i in 1:getncells(grid)
                @test qe.data[i] ≈ 2 * states[i]
            end
            bar(m, u, ∇u, qp_state) = u.u
            a = rand(ndofs(dh2))
            qe_f = QuadPointEvaluator{Float64}(db, bar)
            work!(qe_f, db; a)
            proj = L2Projector(ipu, grid)
            ap = project(proj, qe_f.data, qr)

            @test evaluate_at_grid_nodes(dh2, a, :u) ≈ evaluate_at_grid_nodes(proj, ap)
        end
    end
    @testset "MultiDomain" begin
        qe_state(::QEMat{3}) = rand()
        qe_state(::QEMat{4}) = (rand(), -rand())
        
        ds_left  = DomainSpec(dh1, QEMat{3}(), cvu; set = getcellset(grid, "left"))
        ds_right = DomainSpec(dh1, QEMat{4}(), cvu; set = getcellset(grid, "right"))
        db = setup_domainbuffers(Dict("left" => ds_left, "right" => ds_right))
        states = FerriteAssembly.get_state(db)
        foo(::QEMat{3}, u, ∇u, qp_state) = 3 * qp_state
        foo(::QEMat{4}, u, ∇u, qp_state) = 3 * qp_state[2]
        qe = QuadPointEvaluator{Float64}(db, foo)
        work!(qe, db)
        for (i, s) in states["left"].vals # TODO: Using internals here
            @test qe.data[i] ≈ 3 * s
        end
        for (i, s) in states["right"].vals # TODO: Using internals here
            @test qe.data[i] ≈ 3 * last.(s)
            @test all(first.(s) .≥ 0)
            @test all(last.(s) .≤ 0)
        end

        bar(m, u, ∇u, qp_state) = u
        a = rand(ndofs(dh1))
        qe_f = QuadPointEvaluator{Float64}(db, bar)
        work!(qe_f, db; a)
        proj = L2Projector(ipu, grid)
        ap = project(proj, qe_f.data, qr)

        @test evaluate_at_grid_nodes(dh1, a, :u) ≈ evaluate_at_grid_nodes(proj, ap)
    end
end
