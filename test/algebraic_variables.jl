# Materials and workers used to test algebraic variables (`algebraic_per_cell` in `DomainSpec`)

# Heat equation with the mean value of `u` prescribed to `ū` via the algebraic variable `λ`
struct AVHeatMean{T}
    k::T
    f::T
    ū::T
end
function FerriteAssembly.element_routine!(Ke, re, state, ae, m::AVHeatMean, cv, buffer)
    ru = dof_range(buffer, :u)
    rλ = dof_range(buffer, :λ)
    iλ = only(rλ)
    λ = algebraic_value(FA.get_algebraic_values(buffer, :λ), ae, rλ)
    for q_point in 1:getnquadpoints(cv)
        dΩ = getdetJdV(cv, q_point)
        u = function_value(cv, q_point, ae, ru)
        ∇u = function_gradient(cv, q_point, ae, ru)
        re[iλ] += (u - m.ū) * dΩ
        for (i, I) in pairs(ru)
            δN = shape_value(cv, q_point, i)
            ∇δN = shape_gradient(cv, q_point, i)
            re[I] += (m.k * (∇δN ⋅ ∇u) + δN * (λ - m.f)) * dΩ
            Ke[I, iλ] += δN * dΩ
            Ke[iλ, I] += δN * dΩ
            for (j, J) in pairs(ru)
                Ke[I, J] += m.k * (∇δN ⋅ shape_gradient(cv, q_point, j)) * dΩ
            end
        end
    end
end
function av_heatmean_residual!(re, ae, m, cv, buffer)
    ru = dof_range(buffer, :u)
    rλ = dof_range(buffer, :λ)
    λ = algebraic_value(FA.get_algebraic_values(buffer, :λ), ae, rλ)
    for q_point in 1:getnquadpoints(cv)
        dΩ = getdetJdV(cv, q_point)
        u = function_value(cv, q_point, ae, ru)
        ∇u = function_gradient(cv, q_point, ae, ru)
        re[only(rλ)] += (u - m.ū) * dΩ
        for (i, I) in pairs(ru)
            re[I] += (m.k * (shape_gradient(cv, q_point, i) ⋅ ∇u) + shape_value(cv, q_point, i) * (λ - m.f)) * dΩ
        end
    end
end
FerriteAssembly.element_residual!(re, state, ae, m::AVHeatMean, cv, buffer) = av_heatmean_residual!(re, ae, m, cv, buffer)

# Same as `AVHeatMean`, but only defines the residual, such that automatic differentiation is used
struct AVHeatMeanAD{T}
    m::AVHeatMean{T}
end
FerriteAssembly.element_residual!(re, state, ae, m::AVHeatMeanAD, cv, buffer) = av_heatmean_residual!(re, ae, m.m, cv, buffer)

# Mean value of `u` on the facet domain coupled to the algebraic variable `λ`
struct AVFacetMean{T}
    ū::T
end
function FerriteAssembly.facet_routine!(Ke, re, ae, m::AVFacetMean, fv, buffer)
    ru = dof_range(buffer, :u)
    rλ = dof_range(buffer, :λ)
    iλ = only(rλ)
    λ = algebraic_value(FA.get_algebraic_values(buffer, :λ), ae, rλ)
    for q_point in 1:getnquadpoints(fv)
        dΓ = getdetJdV(fv, q_point)
        u = function_value(fv, q_point, ae, ru)
        re[iλ] += (u - m.ū) * dΓ
        for (i, I) in pairs(ru)
            δN = shape_value(fv, q_point, i)
            re[I] += δN * λ * dΓ
            Ke[I, iλ] += δN * dΓ
            Ke[iλ, I] += δN * dΓ
        end
    end
end
function FerriteAssembly.facet_residual!(re, ae, m::AVFacetMean, fv, buffer)
    ru = dof_range(buffer, :u)
    rλ = dof_range(buffer, :λ)
    λ = algebraic_value(FA.get_algebraic_values(buffer, :λ), ae, rλ)
    for q_point in 1:getnquadpoints(fv)
        dΓ = getdetJdV(fv, q_point)
        re[only(rλ)] += (function_value(fv, q_point, ae, ru) - m.ū) * dΓ
        for (i, I) in pairs(ru)
            re[I] += shape_value(fv, q_point, i) * λ * dΓ
        end
    end
end

# Writes to all entries in `Ke` and `re`, requiring all local dofs to be coupled in the sparsity pattern
struct AVDense end
function FerriteAssembly.element_routine!(Ke, re, state, ae, ::AVDense, cv, buffer)
    fill!(Ke, 1)
    fill!(re, 1)
end
function FerriteAssembly.facet_routine!(Ke, re, ae, ::AVDense, fv, buffer)
    fill!(Ke, 1)
    fill!(re, 1)
end

# Material for testing the inputs to `create_cell_state`, and the caches
struct AVStateMaterial end
FerriteAssembly.create_cell_state(::AVStateMaterial, cv, x, ae, dofrange) = (ae = copy(ae), dofrange = dofrange)
FerriteAssembly.allocate_cell_cache(::AVStateMaterial, cv) = :cellcache
FerriteAssembly.allocate_facet_cache(::AVStateMaterial, fv) = :facetcache

# Record information available in the buffer for each quadrature point (using the `QuadPointEvaluator`)
# (λ, z[1], z[2], global dof for λ, first global dof for z), with `z = (NaN, NaN)` if not available
function av_record(buffer, zbuffer)
    rλ = dof_range(buffer, :λ)
    λ = algebraic_value(FA.get_algebraic_values(buffer, :λ), FA.get_ae(buffer), rλ)
    rz = dof_range(zbuffer, :z)
    z = algebraic_value(FA.get_algebraic_values(zbuffer, :z), FA.get_ae(zbuffer), rz)
    return (λ, z[1], z[2], Float64(celldofs(buffer)[only(rλ)]), Float64(celldofs(zbuffer)[first(rz)]))
end
function FerriteAssembly.eval_quadpoints_cell!(vals::AbstractVector, ::Val{:av_record}, state, ae, m, cv, buffer)
    fill!(vals, av_record(buffer, buffer))
end
function FerriteAssembly.eval_quadpoints_cell!(vals::AbstractVector, ::Val{:av_record_coupled}, state, ae, m, cv, buffer)
    fill!(vals, av_record(buffer, FA.get_coupled_buffer(buffer, :b)))
end

# Worker that can be threaded, but doesn't support domains with dofs shared between all items.
# Counts the number of items worked by the base worker and by the task local workers
struct AVCountWorker
    islocal::Bool
    nbase::Threads.Atomic{Int}
    nlocal::Threads.Atomic{Int}
end
AVCountWorker() = AVCountWorker(false, Threads.Atomic{Int}(0), Threads.Atomic{Int}(0))
FerriteAssembly.can_thread(::AVCountWorker) = true
FerriteAssembly.create_local(w::AVCountWorker) = AVCountWorker(true, w.nbase, w.nlocal)
FerriteAssembly.scatter!(::AVCountWorker, ::AVCountWorker) = nothing
FerriteAssembly.gather!(::AVCountWorker, ::AVCountWorker) = nothing
function FerriteAssembly.work_single_cell!(w::AVCountWorker, cellbuffer)
    Threads.atomic_add!(w.islocal ? w.nlocal : w.nbase, 1)
end
# Same, but supports domains with shared dofs
struct AVSharedCountWorker
    w::AVCountWorker
end
AVSharedCountWorker() = AVSharedCountWorker(AVCountWorker())
FerriteAssembly.can_thread(::AVSharedCountWorker) = true
FerriteAssembly.can_thread_shared_dofs(::AVSharedCountWorker) = true
FerriteAssembly.create_local(w::AVSharedCountWorker) = AVSharedCountWorker(FA.create_local(w.w))
FerriteAssembly.scatter!(::AVSharedCountWorker, ::AVSharedCountWorker) = nothing
FerriteAssembly.gather!(::AVSharedCountWorker, ::AVSharedCountWorker) = nothing
FerriteAssembly.work_single_cell!(w::AVSharedCountWorker, cellbuffer) = FA.work_single_cell!(w.w, cellbuffer)

# Custom cell buffer, that only defines what is required to be worked with the `AVCountWorker`s
mutable struct AVCustomBuffer <: FA.AbstractCellBuffer
    cellid::Int
    Δt::Float64
end
FerriteAssembly.reinit_buffer!(b::AVCustomBuffer, sim, cellnum::Int) = (b.cellid = cellnum)
FerriteAssembly.create_local(b::AVCustomBuffer) = AVCustomBuffer(b.cellid, b.Δt)
FerriteAssembly.get_time_increment(b::AVCustomBuffer) = b.Δt
FerriteAssembly.set_time_increment!(b::AVCustomBuffer, Δt) = (b.Δt = Δt)

# Manual assembly following Ferrite's topic guide on algebraic variables: The local dofs are
# the cell dofs followed by the algebraic dofs.
function av_manual_heatmean!(K, r, dh, cv, m::AVHeatMean, a; set = 1:getncells(dh.grid))
    n = ndofs_per_cell(dh)
    dofs = zeros(Int, n + 1)
    dofs[end] = only(algebraic_dofs(dh, :λ))
    Ke = zeros(n + 1, n + 1)
    re = zeros(n + 1)
    assembler = start_assemble(K, r; fillzero = false)
    for cell in CellIterator(dh, set)
        reinit!(cv, cell)
        copyto!(dofs, 1, celldofs(cell), 1, n)
        ae = a[dofs]
        fill!(Ke, 0)
        fill!(re, 0)
        for q_point in 1:getnquadpoints(cv)
            dΩ = getdetJdV(cv, q_point)
            u = sum(shape_value(cv, q_point, i) * ae[i] for i in 1:n)
            ∇u = sum(shape_gradient(cv, q_point, i) * ae[i] for i in 1:n)
            re[n + 1] += (u - m.ū) * dΩ
            for i in 1:n
                N = shape_value(cv, q_point, i)
                ∇N = shape_gradient(cv, q_point, i)
                re[i] += (m.k * (∇N ⋅ ∇u) - N * m.f + N * ae[n + 1]) * dΩ
                Ke[i, n + 1] += N * dΩ
                Ke[n + 1, i] += N * dΩ
                for j in 1:n
                    Ke[i, j] += m.k * (∇N ⋅ shape_gradient(cv, q_point, j)) * dΩ
                end
            end
        end
        assemble!(assembler, dofs, Ke, re)
    end
    return K, r
end

function av_manual_facetmean!(K, r, dh, fv, m::AVFacetMean, a, facetset)
    n = ndofs_per_cell(dh)
    dofs = zeros(Int, n + 1)
    dofs[end] = only(algebraic_dofs(dh, :λ))
    Ke = zeros(n + 1, n + 1)
    re = zeros(n + 1)
    assembler = start_assemble(K, r; fillzero = false)
    for facet in FacetIterator(dh, facetset)
        reinit!(fv, facet)
        copyto!(dofs, 1, celldofs(facet), 1, n)
        ae = a[dofs]
        fill!(Ke, 0)
        fill!(re, 0)
        for q_point in 1:getnquadpoints(fv)
            dΓ = getdetJdV(fv, q_point)
            u = sum(shape_value(fv, q_point, i) * ae[i] for i in 1:n)
            re[n + 1] += (u - m.ū) * dΓ
            for i in 1:n
                N = shape_value(fv, q_point, i)
                re[i] += N * ae[n + 1] * dΓ
                Ke[i, n + 1] += N * dΓ
                Ke[n + 1, i] += N * dΓ
            end
        end
        assemble!(assembler, dofs, Ke, re)
    end
    return K, r
end

# Reverse the dof numbering, such that the algebraic dofs are numbered first. Returns the permuted `a`
function av_renumber!(dh, a)
    perm = collect(ndofs(dh):-1:1)
    renumber!(dh, perm)
    a_new = similar(a)
    a_new[perm] .= a
    return a_new
end

# Work a single item, used for testing that this doesn't allocate
function av_work_item!(worker, itembuffer, sim, item)
    FA.reinit_buffer!(itembuffer, sim, item)
    FA.work_single!(worker, itembuffer)
    return nothing
end

@testset "Algebraic variables" begin
    ip = Lagrange{RefQuadrilateral, 1}()
    cv = CellValues(QuadratureRule{RefQuadrilateral}(2), ip)
    fv = FacetValues(FacetQuadratureRule{RefQuadrilateral}(2), ip)
    function setup_dh(nel = (6, 5); variables = (:λ => AlgebraicVariable(),))
        grid = generate_grid(Quadrilateral, nel)
        dh = DofHandler(grid)
        add!(dh, :u, ip)
        foreach(v -> add!(dh, v.first, v.second), variables)
        close!(dh)
        return dh
    end
    material = AVHeatMean(1.2, 0.7, 2.5)
    facetmaterial = AVFacetMean(0.4)

    @testset "Layout computation" begin
        # Only standalone variables: The cell dofs are not moved
        r = FA.compute_fieldwise_ranges([:u => 8, :p => 4], Dict{Symbol, Vector{Pair{Symbol, Int}}}(), [:a => 1, :b => 3])
        @test r.field_sources == [1:8, 9:12]
        @test r.field_destinations == [1, 9]
        @test r.field_ranges == [1:8, 9:12]
        @test r.variable_names == [:a, :b]
        @test r.variable_ranges == [13:13, 14:16]
        @test r.ndofs == 16
        # Variables attached to a field are located directly after the field's cell dofs
        r = FA.compute_fieldwise_ranges([:u => 8, :p => 4], Dict(:u => [:H => 3, :γ => 2]), [:a => 1])
        @test r.field_sources == [1:8, 9:12]
        @test r.field_destinations == [1, 14]
        @test r.field_ranges == [1:13, 14:17]
        @test r.variable_names == [:H, :γ, :a]
        @test r.variable_ranges == [9:11, 12:13, 18:18]
        @test r.ndofs == 18
        r = FA.compute_fieldwise_ranges([:u => 8, :p => 4], Dict(:p => [:H => 2]), Pair{Symbol, Int}[])
        @test r.field_destinations == [1, 9]
        @test r.field_ranges == [1:8, 9:14]
        @test r.variable_ranges == [13:14]
        @test r.ndofs == 14
    end

    @testset "Ferrite helpers" begin
        dh = setup_dh(; variables = (:λ => AlgebraicVariable(), :z => AlgebraicVariable{Vec{2}}()))
        @test FA._algebraic_names(dh) == [:λ, :z]
        @test FA._algebraic_names(dh.subdofhandlers[1]) == [:λ, :z]
        @test isempty(FA._algebraic_names(setup_dh(; variables = ())))
        @test FA._algebraic_variable(dh, :λ) isa AlgebraicVariable
        @test getnbasefunctions(AlgebraicValues(FA._algebraic_variable(dh, :z))) == 2
        @test_throws ArgumentError FA._algebraic_variable(dh, :u)
    end

    @testset "Setup and accessors" begin
        dh = setup_dh(; variables = (:λ => AlgebraicVariable(), :z => AlgebraicVariable{Vec{2}}()))
        a = rand(ndofs(dh))
        n = ndofs_per_cell(dh)
        for (db, set) in (
                (setup_domainbuffer(DomainSpec(dh, nothing, cv; algebraic_per_cell = (:z, :λ))), 1:getncells(dh.grid)),
                (setup_domainbuffer(DomainSpec(dh, nothing, fv; algebraic_per_cell = (:z, :λ), set = getfacetset(dh.grid, "left"))), getfacetset(dh.grid, "left")),
            )
            buffer = FA.get_itembuffer(db)
            @test Ferrite.getfieldnames(buffer) == (:u, :z, :λ)
            @test dof_range(buffer, :u) == 1:n
            @test dof_range(buffer, :z) == (n + 1):(n + 2)
            @test dof_range(buffer, :λ) == (n + 3):(n + 3)
            @test FA.get_algebraic_values(buffer, :z) isa AlgebraicValues
            @test_throws Exception FA.get_algebraic_values(buffer, :u)
            @test size(FA.get_Ke(buffer)) == (n + 3, n + 3)
            @test length(FA.get_re(buffer)) == length(FA.get_ae(buffer)) == length(FA.get_aeold(buffer)) == n + 3
            sim = Simulation(db, a)
            item = first(FA.get_set(db))
            FA.reinit_buffer!(buffer, sim, item)
            @test celldofs(buffer) == vcat(celldofs(dh, first(item)), algebraic_dofs(dh, :z), algebraic_dofs(dh, :λ))
            @test FA.get_ae(buffer) == a[celldofs(buffer)]
            @test all(isnan, FA.get_aeold(buffer))
            z = algebraic_value(FA.get_algebraic_values(buffer, :z), FA.get_ae(buffer), dof_range(buffer, :z))
            @test z == algebraic_value(dh, a, :z)
            λ = algebraic_value(FA.get_algebraic_values(buffer, :λ), FA.get_ae(buffer), dof_range(buffer, :λ))
            @test λ == algebraic_value(dh, a, :λ)
        end
        # No algebraic variables
        buffer = FA.get_itembuffer(setup_domainbuffer(DomainSpec(dh, nothing, cv)))
        @test Ferrite.getfieldnames(buffer) == (:u,)
        @test FA.get_layout(buffer) === nothing
        @test_throws ArgumentError FA.get_algebraic_values(buffer, :λ)
        @test length(celldofs(buffer)) == n
    end

    @testset "Validation" begin
        dh = setup_dh()
        @test_throws ArgumentError DomainSpec(dh, nothing, cv; algebraic_per_cell = :λ)
        @test_throws ArgumentError DomainSpec(dh, nothing, cv; algebraic_per_cell = ("λ",))
        @test_throws ArgumentError setup_domainbuffer(DomainSpec(dh, nothing, cv; algebraic_per_cell = (:u,)))
        @test_throws ArgumentError setup_domainbuffer(DomainSpec(dh, nothing, cv; algebraic_per_cell = (:nonexistent,)))
        @test_throws ArgumentError setup_domainbuffer(DomainSpec(dh, nothing, cv; algebraic_per_cell = (:λ, :λ)))
        @test_throws ArgumentError setup_domainbuffer(DomainSpec(dh, nothing, fv; algebraic_per_cell = (:nonexistent,), set = getfacetset(dh.grid, "left")))
        @test_throws ArgumentError setup_domainbuffers(Dict("a" => DomainSpec(dh, nothing, cv; algebraic_per_cell = (:u,))))
    end

    @testset "Cell domain" begin
        dh = setup_dh()
        a = rand(ndofs(dh))
        db = setup_domainbuffer(DomainSpec(dh, material, cv; algebraic_per_cell = (:λ,)))
        couplings = FA.algebraic_couplings(db)
        @test couplings isa Vector{Union{CellCoupling, FacetCoupling}}
        @test only(couplings) isa CellCoupling
        K_ref = allocate_matrix(dh; algebraic_couplings = couplings)
        r_ref = zeros(ndofs(dh))
        av_manual_heatmean!(K_ref, r_ref, dh, cv, material, a)

        K = allocate_matrix(dh; algebraic_couplings = couplings)
        r = zeros(ndofs(dh))
        work!(KeReAssembler(K, r), db; a = a)
        @test K ≈ K_ref
        @test r ≈ r_ref

        K2 = allocate_matrix(dh; algebraic_couplings = couplings)
        r2 = zeros(ndofs(dh))
        work!(start_assemble(K2, r2), db; a = a)
        @test K2 ≈ K_ref
        @test r2 ≈ r_ref

        r3 = zeros(ndofs(dh))
        work!(ReAssembler(r3), Simulation(db, a))
        @test r3 ≈ r_ref

        # Solve the (linear) problem and check that the mean value is as prescribed
        a .-= K \ r
        work!(ReAssembler(r3), db; a = a)
        @test norm(r3) < 1.0e-10
        integrator = SimpleIntegrator((u, ∇u, state) -> (1.0, u), (0.0, 0.0))
        work!(integrator, db; a = a)
        @test integrator.val[2] / integrator.val[1] ≈ material.ū

        # Automatic differentiation
        for autodiffbuffer in (false, true)
            db_ad = setup_domainbuffer(DomainSpec(dh, AVHeatMeanAD(material), cv; algebraic_per_cell = (:λ,)); autodiffbuffer)
            @test (FA.get_itembuffer(db_ad) isa FA.AutoDiffCellBuffer) == autodiffbuffer
            K_ad = allocate_matrix(dh; algebraic_couplings = FA.algebraic_couplings(db_ad))
            r_ad = zeros(ndofs(dh))
            a_ad = rand(ndofs(dh))
            work!(KeReAssembler(K_ad, r_ad), db_ad; a = a_ad)
            @test K_ad ≈ K_ref
            fill!(r_ref, 0)
            work!(ReAssembler(r_ref), db; a = a_ad)
            @test r_ad ≈ r_ref
        end
    end

    @testset "Facet domain" begin
        dh = setup_dh()
        a = rand(ndofs(dh))
        facetset = getfacetset(dh.grid, "right")
        db = setup_domainbuffer(DomainSpec(dh, facetmaterial, fv; set = facetset, algebraic_per_cell = (:λ,)))
        couplings = FA.algebraic_couplings(db)
        @test only(couplings) isa FacetCoupling
        K_ref = allocate_matrix(dh; algebraic_couplings = couplings)
        r_ref = zeros(ndofs(dh))
        av_manual_facetmean!(K_ref, r_ref, dh, fv, facetmaterial, a, facetset)
        @test norm(r_ref) > 0.1

        K = allocate_matrix(dh; algebraic_couplings = couplings)
        r = zeros(ndofs(dh))
        work!(KeReAssembler(K, r), db; a = a)
        @test K ≈ K_ref
        @test r ≈ r_ref

        r2 = zeros(ndofs(dh))
        work!(ReAssembler(r2), db; a = a)
        @test r2 ≈ r_ref
    end

    @testset "Variable in both cell and facet domains" begin
        dh = setup_dh()
        a = rand(ndofs(dh))
        facetset = getfacetset(dh.grid, "top")
        cellbuffers = setup_domainbuffers(Dict(
            "left" => DomainSpec(dh, material, cv; algebraic_per_cell = (:λ,), set = 1:10),
            "right" => DomainSpec(dh, material, cv; algebraic_per_cell = (:λ,), set = 11:getncells(dh.grid)));
            )
        facetbuffer = setup_domainbuffer(DomainSpec(dh, facetmaterial, fv; set = facetset, algebraic_per_cell = (:λ,)))
        couplings = FA.algebraic_couplings(cellbuffers, Simulation(facetbuffer))
        @test length(couplings) == 3
        @test count(c -> c isa CellCoupling, couplings) == 2
        @test count(c -> c isa FacetCoupling, couplings) == 1

        K_ref = allocate_matrix(dh; algebraic_couplings = couplings)
        r_ref = zeros(ndofs(dh))
        av_manual_heatmean!(K_ref, r_ref, dh, cv, material, a)
        av_manual_facetmean!(K_ref, r_ref, dh, fv, facetmaterial, a, facetset)

        K = allocate_matrix(dh; algebraic_couplings = couplings)
        r = zeros(ndofs(dh))
        work!(KeReAssembler(K, r), cellbuffers; a = a)
        work!(KeReAssembler(K, r; fillzero = false), facetbuffer; a = a)
        @test K ≈ K_ref
        @test r ≈ r_ref
    end

    @testset "Sparsity pattern" begin
        # Multiple fields and algebraic variables: All local dofs must be coupled
        grid = generate_grid(Quadrilateral, (3, 3))
        dh = DofHandler(grid)
        add!(dh, :u, ip^2)
        add!(dh, :p, ip)
        add!(dh, :z, AlgebraicVariable{Vec{2}}())
        add!(dh, :λ, AlgebraicVariable())
        add!(dh, :unused, AlgebraicVariable())
        close!(dh)
        cvs = (u = CellValues(QuadratureRule{RefQuadrilateral}(2), ip^2), p = cv)
        fvs = (u = FacetValues(FacetQuadratureRule{RefQuadrilateral}(2), ip^2), p = fv)
        cellbuffer = setup_domainbuffer(DomainSpec(dh, AVDense(), cvs; algebraic_per_cell = (:λ, :z)))
        facetbuffer = setup_domainbuffer(DomainSpec(dh, AVDense(), fvs; algebraic_per_cell = (:z,), set = getfacetset(grid, "left")))
        nodomain = setup_domainbuffer(DomainSpec(dh, AVDense(), cvs))
        emptydomain = setup_domainbuffer(DomainSpec(dh, AVDense(), cvs; algebraic_per_cell = (:λ,), set = Int[]))
        @test isempty(FA.algebraic_couplings(nodomain, emptydomain))
        @test isempty(FA.algebraic_couplings())

        ncells = getncells(grid)
        nfacets = length(getfacetset(grid, "left"))
        λdof = only(algebraic_dofs(dh, :λ))
        zdofs = algebraic_dofs(dh, :z)
        # Assembly throws if entries are missing in the sparsity pattern
        K = allocate_matrix(dh; algebraic_couplings = FA.algebraic_couplings(cellbuffer))
        r = zeros(ndofs(dh))
        work!(KeReAssembler(K, r), cellbuffer)
        @test K[λdof, λdof] == ncells
        @test K[λdof, zdofs[2]] == K[zdofs[1], λdof] == K[zdofs[1], zdofs[2]] == ncells
        @test r[λdof] == r[zdofs[1]] == ncells
        @test all(!iszero, K[λdof, 1:(first(zdofs) - 1)]) # Coupled to all fields in all cells
        K = allocate_matrix(dh; algebraic_couplings = FA.algebraic_couplings(cellbuffer, facetbuffer))
        work!(KeReAssembler(K, r), cellbuffer)
        work!(KeReAssembler(K, r; fillzero = false), facetbuffer)
        @test K[zdofs[1], zdofs[2]] == ncells + nfacets
        @test K[λdof, zdofs[2]] == ncells
        @test r[zdofs[1]] == ncells + nfacets
    end

    @testset "Constraints" begin
        # λ is prescribed, and one dof for u is constrained to the algebraic variable H, which
        # is not used in any domain
        dh = setup_dh(; variables = (:λ => AlgebraicVariable(), :H => AlgebraicVariable()))
        λdof = only(algebraic_dofs(dh, :λ))
        Hdof = only(algebraic_dofs(dh, :H))
        udof = celldofs(dh, getncells(dh.grid))[3]
        ch = ConstraintHandler(dh)
        add!(ch, Dirichlet(:u, getfacetset(dh.grid, "left"), Returns(0.3)))
        add!(ch, AffineConstraint(λdof, Pair{Int, Float64}[], -0.8))
        add!(ch, AffineConstraint(udof, [Hdof => 2.0], 0.1))
        close!(ch)
        db = setup_domainbuffer(DomainSpec(dh, material, cv; algebraic_per_cell = (:λ,)))
        K = allocate_matrix(dh, ch; algebraic_couplings = FA.algebraic_couplings(db))
        r = zeros(ndofs(dh))

        # Local application of constraints
        a_local = zeros(ndofs(dh))
        apply!(a_local, ch)
        work!(KeReAssembler(K, r; ch, apply_zero = true), db; a = a_local)
        a_local .-= K \ r
        apply!(a_local, ch)

        # Global application of constraints
        a_global = zeros(ndofs(dh))
        apply!(a_global, ch)
        work!(KeReAssembler(K, r), db; a = a_global)
        apply_zero!(K, r, ch)
        a_global .-= K \ r
        apply!(a_global, ch)

        @test a_local ≈ a_global
        @test a_local[λdof] ≈ -0.8
        @test a_local[udof] ≈ 2 * a_local[Hdof] + 0.1
        @test abs(a_local[Hdof]) > 1.0e-3
        # The solution fulfills the equations for the free dofs
        work!(ReAssembler(r), db; a = a_local)
        @test norm(r[ch.free_dofs]) < 1.0e-10
        @test abs(r[udof]) < 1.0e-10 # The condensed equation for H is 2*r[udof] = 0
        @test abs(r[λdof]) > 1.0e-3  # λ is prescribed
    end

    @testset "Workers and scaling" begin
        dh = setup_dh()
        dh_plain = setup_dh(; variables = ())
        a = rand(ndofs(dh))
        a_plain = a[1:ndofs(dh_plain)] # Algebraic dofs are numbered last
        fun = (u, ∇u, state) -> (u, ∇u)
        val0 = (0.0, zero(Vec{2}))
        for (values, kwargs) in ((cv, (;)), (fv, (; set = getfacetset(dh.grid, "left"))))
            db = setup_domainbuffer(DomainSpec(dh, material, values; algebraic_per_cell = (:λ,), kwargs...))
            db_plain = setup_domainbuffer(DomainSpec(dh_plain, material, values; kwargs...))
            integrator = SimpleIntegrator(fun, val0)
            work!(integrator, db; a = a)
            integrator_plain = SimpleIntegrator(fun, val0)
            work!(integrator_plain, db_plain; a = a_plain)
            @test integrator.val[1] ≈ integrator_plain.val[1]
            @test integrator.val[2] ≈ integrator_plain.val[2]
        end
        db = setup_domainbuffer(DomainSpec(dh, material, cv; algebraic_per_cell = (:λ,)))
        db_plain = setup_domainbuffer(DomainSpec(dh_plain, material, cv))
        qe = QuadPointEvaluator{Float64}(db, (m, u, ∇u, state) -> u + ∇u[1])
        work!(qe, db; a = a)
        qe_plain = QuadPointEvaluator{Float64}(db_plain, (m, u, ∇u, state) -> u + ∇u[1])
        work!(qe_plain, db_plain; a = a_plain)
        @test qe.data == qe_plain.data

        # Scaling
        @test Set(keys(ElementResidualScaling(dh).factors)) == Set((:u, :λ))
        @test Set(keys(ElementResidualScaling(dh.subdofhandlers[1]).factors)) == Set((:u, :λ))
        @test Set(keys(ElementResidualScaling(dh_plain).factors)) == Set((:u,))
        @test Set(keys(ElementResidualScaling(dh_plain.subdofhandlers[1]).factors)) == Set((:u,))
        # Reference values
        factor_u = 0.0
        factor_λ = 0.0
        buffer = FA.get_itembuffer(db)
        sim = Simulation(db, a)
        for cellnr in 1:getncells(dh.grid)
            FA.reinit_buffer!(buffer, sim, cellnr)
            re = FA.get_re(buffer)
            FA.element_residual!(re, nothing, FA.get_ae(buffer), material, cv, buffer)
            factor_u += norm(re[dof_range(buffer, :u)])
            factor_λ += norm(re[dof_range(buffer, :λ)])
        end
        r = zeros(ndofs(dh))
        scaling = ElementResidualScaling(dh)
        work!(ReAssembler(r; scaling), db; a = a)
        @test scaling.factors[:u] ≈ factor_u
        @test scaling.factors[:λ] ≈ factor_λ
        for autodiffbuffer in (false, true)
            db_ad = setup_domainbuffer(DomainSpec(dh, AVHeatMeanAD(material), cv; algebraic_per_cell = (:λ,)); autodiffbuffer)
            K = allocate_matrix(dh; algebraic_couplings = FA.algebraic_couplings(db_ad))
            scaling = ElementResidualScaling(dh, 1)
            work!(KeReAssembler(K, r; scaling), db_ad; a = a)
            @test scaling.factors[:λ] ≈ factor_λ
        end
    end

    @testset "States and caches" begin
        dh = setup_dh(; variables = (:λ => AlgebraicVariable(), :z => AlgebraicVariable{Vec{2}}()))
        a = rand(ndofs(dh))
        n = ndofs_per_cell(dh)
        db = setup_domainbuffer(DomainSpec(dh, AVStateMaterial(), cv; algebraic_per_cell = (:z, :λ)); a = a)
        for cellnr in (1, getncells(dh.grid))
            state = FA.get_state(db, cellnr)
            @test keys(state.dofrange) == (:u, :z, :λ)
            @test state.dofrange == (u = 1:n, z = (n + 1):(n + 2), λ = (n + 3):(n + 3))
            @test state.ae == a[vcat(celldofs(dh, cellnr), algebraic_dofs(dh, :z), algebraic_dofs(dh, :λ))]
            @test FA.get_old_state(db, cellnr).ae == state.ae
        end
        @test FA.get_user_cache(FA.get_itembuffer(db)) === :cellcache
        db = setup_domainbuffer(DomainSpec(dh, AVStateMaterial(), cv; algebraic_per_cell = (:z, :λ)))
        @test all(isnan, FA.get_state(db, 1).ae)
        @test length(FA.get_state(db, 1).ae) == n + 3
        db = setup_domainbuffer(DomainSpec(dh, AVStateMaterial(), fv; algebraic_per_cell = (:z, :λ), set = getfacetset(dh.grid, "left")))
        @test FA.get_user_cache(FA.get_itembuffer(db)) === :facetcache
    end

    @testset "Renumbering, material replacement, and empty domains" begin
        variables = (:λ => AlgebraicVariable(), :z => AlgebraicVariable{Vec{2}}())
        dh = setup_dh(; variables)
        a = rand(ndofs(dh))
        # Setup all buffers before renumbering
        db = setup_domainbuffer(DomainSpec(dh, material, cv; algebraic_per_cell = (:λ,)))
        db_record = setup_domainbuffer(DomainSpec(dh, material, cv; algebraic_per_cell = (:λ, :z)))
        db_facet = setup_domainbuffer(DomainSpec(dh, facetmaterial, fv; algebraic_per_cell = (:λ,), set = getfacetset(dh.grid, "left")))
        threaded = [setup_domainbuffer(DomainSpec(dh, material, cv; algebraic_per_cell = (:λ, :z)); threading = true, autodiffbuffer)
            for autodiffbuffer in (false, true)]
        qe = QuadPointEvaluator{NTuple{5, Float64}}(db_record, :av_record)
        work!(qe, db_record; a = a)
        λ0 = algebraic_value(dh, a, :λ)
        z0 = algebraic_value(dh, a, :z)
        @test all(v -> all(==((λ0, z0[1], z0[2], ndofs(dh) - 2.0, ndofs(dh) - 1.0)), v), qe.data)

        a = av_renumber!(dh, a)
        @test algebraic_dofs(dh, :λ) == [3]
        @test algebraic_dofs(dh, :z) == [2, 1]
        @test algebraic_value(dh, a, :z) == z0
        expected = (λ0, z0[1], z0[2], 3.0, 2.0)
        work!(qe, db_record; a = a)
        @test all(v -> all(==(expected), v), qe.data)
        for db_threaded in threaded
            qe_threaded = QuadPointEvaluator{NTuple{5, Float64}}(db_threaded, :av_record)
            work!(qe_threaded, db_threaded; a = a)
            @test all(v -> all(==(expected), v), qe_threaded.data)
            foreach(FA.get_locals(FA.get_itembuffer(db_threaded))) do buffer
                @test celldofs(buffer)[dof_range(buffer, :z)] == [2, 1]
            end
        end

        # Assembly with buffers created before and after renumbering
        for (old, new) in (
                (db, setup_domainbuffer(DomainSpec(dh, material, cv; algebraic_per_cell = (:λ,)))),
                (db_facet, setup_domainbuffer(DomainSpec(dh, facetmaterial, fv; algebraic_per_cell = (:λ,), set = getfacetset(dh.grid, "left")))),
            )
            K_new = allocate_matrix(dh; algebraic_couplings = FA.algebraic_couplings(new))
            r_new = zeros(ndofs(dh))
            work!(KeReAssembler(K_new, r_new), new; a = a)
            K_old = allocate_matrix(dh; algebraic_couplings = FA.algebraic_couplings(old))
            r_old = zeros(ndofs(dh))
            work!(KeReAssembler(K_old, r_old), old; a = a)
            @test K_old == K_new
            @test r_old == r_new
            @test norm(r_old) > 0.1
        end

        # Replacing the material
        material2 = AVHeatMean(2 * material.k, material.f, material.ū - 1)
        db2 = FA.replace_material(db, Returns(material2))
        db2_ref = setup_domainbuffer(DomainSpec(dh, material2, cv; algebraic_per_cell = (:λ,)))
        r2 = zeros(ndofs(dh))
        r2_ref = zeros(ndofs(dh))
        work!(ReAssembler(r2), db2; a = a)
        work!(ReAssembler(r2_ref), db2_ref; a = a)
        @test r2 == r2_ref

        # Empty domains
        for threading in (false, true)
            db_empty = setup_domainbuffer(DomainSpec(dh, material, cv; algebraic_per_cell = (:λ,), set = Int[]); threading)
            worker = AVSharedCountWorker()
            work!(worker, db_empty; a = a)
            @test worker.w.nbase[] == worker.w.nlocal[] == 0
            @test isempty(FA.algebraic_couplings(db_empty))
        end
    end

    @testset "Coupled simulations" begin
        grid = generate_grid(Quadrilateral, (4, 4))
        dh_a = DofHandler(grid)
        add!(dh_a, :u, ip)
        add!(dh_a, :λ, AlgebraicVariable())
        close!(dh_a)
        dh_b = DofHandler(grid)
        add!(dh_b, :u, ip)
        add!(dh_b, :other, AlgebraicVariable())
        add!(dh_b, :z, AlgebraicVariable{Vec{2}}())
        close!(dh_b)
        a_a = rand(ndofs(dh_a))
        a_b = rand(ndofs(dh_b))
        λ0 = algebraic_value(dh_a, a_a, :λ)
        z0 = algebraic_value(dh_b, a_b, :z)
        for threading in (false, true), ad_a in (false, true), ad_b in (false, true)
            sim_a = Simulation(setup_domainbuffer(DomainSpec(dh_a, material, cv; algebraic_per_cell = (:λ,)); threading, autodiffbuffer = ad_a), a_a)
            sim_b = Simulation(setup_domainbuffer(DomainSpec(dh_b, material, cv; algebraic_per_cell = (:z,)); threading, autodiffbuffer = ad_b), a_b)
            group = CoupledSimulations((a = sim_a,); refs = (b = sim_b,))
            qe = QuadPointEvaluator{NTuple{5, Float64}}(group.a.db, :av_record_coupled)
            work!(qe, group.a)
            expected = (λ0, z0[1], z0[2], Float64(only(algebraic_dofs(dh_a, :λ))), Float64(first(algebraic_dofs(dh_b, :z))))
            @test all(v -> all(==(expected), v), qe.data)
        end
        # Renumbering the partner's dofhandler after setting up the simulations
        sims = map(((false, false), (true, false), (true, true))) do (threading, autodiffbuffer)
            sim_a = Simulation(setup_domainbuffer(DomainSpec(dh_a, material, cv; algebraic_per_cell = (:λ,)); threading, autodiffbuffer), a_a)
            db_b = setup_domainbuffer(DomainSpec(dh_b, material, cv; algebraic_per_cell = (:z,)); threading, autodiffbuffer)
            return (sim_a, db_b)
        end
        a_b = av_renumber!(dh_b, a_b)
        @test algebraic_dofs(dh_b, :z) == [2, 1]
        expected = (λ0, z0[1], z0[2], Float64(only(algebraic_dofs(dh_a, :λ))), 2.0)
        for (sim_a, db_b) in sims
            group = CoupledSimulations((a = sim_a,); refs = (b = Simulation(db_b, a_b),))
            qe = QuadPointEvaluator{NTuple{5, Float64}}(group.a.db, :av_record_coupled)
            work!(qe, group.a)
            @test all(v -> all(==(expected), v), qe.data)
        end
    end

    @testset "Threading" begin
        dh = setup_dh((20, 20))
        a = rand(ndofs(dh))
        ncells = getncells(dh.grid)
        facetset = getfacetset(dh.grid, "right")
        db = setup_domainbuffer(DomainSpec(dh, material, cv; algebraic_per_cell = (:λ,)))
        dbt = setup_domainbuffer(DomainSpec(dh, material, cv; algebraic_per_cell = (:λ,)); threading = true)
        dbf = setup_domainbuffer(DomainSpec(dh, facetmaterial, fv; algebraic_per_cell = (:λ,), set = facetset))
        dbft = setup_domainbuffer(DomainSpec(dh, facetmaterial, fv; algebraic_per_cell = (:λ,), set = facetset); threading = true)
        couplings = FA.algebraic_couplings(db, dbf)

        @testset "traits" begin
            K = allocate_matrix(dh; algebraic_couplings = couplings)
            r = zeros(ndofs(dh))
            @test FA.can_thread_shared_dofs(KeReAssembler(K, r; atomic = true))
            @test !FA.can_thread_shared_dofs(KeReAssembler(K, r))
            @test !FA.can_thread_shared_dofs(KeReAssembler(start_assemble(K, r; atomic = true)))
            @test !FA.can_thread_shared_dofs(start_assemble(K, r; atomic = true))
            @test !FA.can_thread_shared_dofs(ReAssembler(r))
            @test FA.can_thread_shared_dofs(SimpleIntegrator(Returns(1.0), 0.0))
            @test FA.can_thread_shared_dofs(QuadPointEvaluator{Float64}(db, Returns(1.0)))
            @test !FA.can_thread_shared_dofs(Integrator(zeros(1)))
            @test !FA.can_thread_shared_dofs(AVCountWorker())
            # The task-local assemblers keep the atomic mode
            assembler = KeReAssembler(K, r; atomic = true)
            @test FA.create_local(assembler) isa typeof(assembler)
            @test typeof(FA.create_local(assembler).a) === typeof(start_assemble(K, r; atomic = true))
            @test typeof(FA.create_local(KeReAssembler(K, r)).a) === typeof(start_assemble(K, r))
        end

        @testset "task local buffers" begin
            for autodiffbuffer in (false, true)
                dbt_ad = setup_domainbuffer(DomainSpec(dh, AVHeatMeanAD(material), cv; algebraic_per_cell = (:λ,)); threading = true, autodiffbuffer)
                itembuffers = FA.get_itembuffer(dbt_ad)
                base = FA.get_base(itembuffers)
                locals = FA.get_locals(itembuffers)
                for buffer in locals
                    @test FA.get_layout(buffer) === FA.get_layout(base)
                    @test Ferrite.getfieldnames(buffer) == (:u, :λ)
                    @test celldofs(buffer) !== celldofs(base)
                    @test celldofs(buffer)[dof_range(buffer, :λ)] == algebraic_dofs(dh, :λ)
                end
                @test allunique(objectid.(celldofs.(locals)))
            end
        end

        @testset "assembly" begin
            for (seq, thr) in ((db, dbt), (dbf, dbft))
                K_ref = allocate_matrix(dh; algebraic_couplings = couplings)
                r_ref = zeros(ndofs(dh))
                work!(KeReAssembler(K_ref, r_ref), seq; a = a)

                # Atomic assembly is done multithreaded
                K = allocate_matrix(dh; algebraic_couplings = couplings)
                r = zeros(ndofs(dh))
                @test_logs work!(KeReAssembler(K, r; atomic = true), thr; a = a)
                @test K ≈ K_ref
                @test r ≈ r_ref

                # Workers that are not safe for shared dofs fall back to sequential assembly
                @test_logs (:warn, r"algebraic variables") work!(KeReAssembler(K, r), thr; a = a)
                @test K ≈ K_ref
                @test r ≈ r_ref
                fill!(r, 0)
                @test_logs (:warn, r"algebraic variables") work!(ReAssembler(r), thr; a = a)
                @test r ≈ r_ref
                @test_logs (:warn, r"algebraic variables") work!(start_assemble(K, r), thr; a = a)
                @test K ≈ K_ref
                @test r ≈ r_ref
            end
        end

        @testset "choice of sequential or threaded work per domain" begin
            specs = Dict(
                "with" => DomainSpec(dh, material, cv; algebraic_per_cell = (:λ,), set = 1:100),
                "without" => DomainSpec(dh, material, cv; set = 101:ncells))
            buffers = setup_domainbuffers(specs; threading = true)
            worker = AVCountWorker()
            @test_logs (:warn, r"AVCountWorker") work!(worker, buffers)
            @test worker.nbase[] == 100
            @test worker.nlocal[] == ncells - 100
            worker = AVSharedCountWorker()
            @test_logs work!(worker, buffers)
            @test worker.w.nbase[] == 0
            @test worker.w.nlocal[] == ncells
            worker = AVCountWorker()
            @test_logs (:warn, r"AVCountWorker") work!(worker, buffers["with"])
            @test (worker.nbase[], worker.nlocal[]) == (100, 0)
            worker = AVCountWorker()
            @test_logs work!(worker, buffers["without"])
            @test (worker.nbase[], worker.nlocal[]) == (0, ncells - 100)
            worker = AVSharedCountWorker()
            @test_logs work!(worker, buffers["with"])
            @test (worker.w.nbase[], worker.w.nlocal[]) == (0, 100)

            # Integration still works threaded and gives the same result
            integrator = SimpleIntegrator((u, ∇u, state) -> u, 0.0)
            work!(integrator, db; a = a)
            integrator_threaded = SimpleIntegrator((u, ∇u, state) -> u, 0.0)
            @test_logs work!(integrator_threaded, dbt; a = a)
            @test integrator_threaded.val ≈ integrator.val
        end
    end

    @testset "Custom item buffers without algebraic variables" begin
        dh = setup_dh(; variables = ())
        ncells = getncells(dh.grid)
        set = collect(1:ncells)
        states = FA.StateVariables(Dict{Int, Nothing}(), Dict{Int, Nothing}())
        sdh = dh.subdofhandlers[1]
        db = FA.DomainBuffer(set, AVCustomBuffer(0, 0.0), states, sdh)
        @test !FA.has_algebraic_dofs(FA.get_itembuffer(db))
        worker = AVCountWorker()
        work!(worker, db)
        @test (worker.nbase[], worker.nlocal[]) == (ncells, 0)
        dbt = FA.ThreadedDomainBuffer(set, AVCustomBuffer(0, 0.0), states, sdh)
        worker = AVCountWorker()
        @test_logs work!(worker, dbt)
        @test (worker.nbase[], worker.nlocal[]) == (0, ncells)
        @test isempty(FA.algebraic_couplings(db, dbt))
    end

    @testset "Allocations" begin
        function setup_assembly(nel)
            dh = setup_dh(nel)
            db = setup_domainbuffer(DomainSpec(dh, material, cv; algebraic_per_cell = (:λ,)))
            dbf = setup_domainbuffer(DomainSpec(dh, facetmaterial, fv; algebraic_per_cell = (:λ,), set = getfacetset(dh.grid, "left")))
            K = allocate_matrix(dh; algebraic_couplings = FA.algebraic_couplings(db, dbf))
            return db, dbf, KeReAssembler(K, zeros(ndofs(dh))), zeros(ndofs(dh))
        end
        work_allocations(worker, db, a) = @allocated work!(worker, Simulation(db, a))
        item_allocations(worker, buffer, sim, item) = @allocated av_work_item!(worker, buffer, sim, item)

        for i in (1, 2) # cell and facet
            domain_small, worker_small, a_small = setup_assembly((4, 4))[[i, 3, 4]]
            domain_large, worker_large, a_large = setup_assembly((16, 16))[[i, 3, 4]]
            # Working an item doesn't allocate
            sim = Simulation(domain_small, a_small)
            buffer = FA.get_itembuffer(domain_small)
            item = first(FA.get_set(domain_small))
            item_allocations(worker_small, buffer, sim, item)
            @test item_allocations(worker_small, buffer, sim, item) == 0
            # Updating the algebraic dofs allocates once before working a domain, independent of the number of items
            work_allocations(worker_small, domain_small, a_small)
            work_allocations(worker_large, domain_large, a_large)
            @test work_allocations(worker_small, domain_small, a_small) == work_allocations(worker_large, domain_large, a_large)
        end

        # Task-local buffers
        dh = setup_dh((4, 4))
        a = zeros(ndofs(dh))
        dbt = setup_domainbuffer(DomainSpec(dh, material, cv; algebraic_per_cell = (:λ,)); threading = true)
        K = allocate_matrix(dh; algebraic_couplings = FA.algebraic_couplings(dbt))
        worker = FA.create_local(KeReAssembler(K, zeros(ndofs(dh)); atomic = true))
        sim = Simulation(dbt, a)
        buffer = first(FA.get_locals(FA.get_itembuffer(dbt)))
        item_allocations(worker, buffer, sim, 1)
        @test item_allocations(worker, buffer, sim, 1) == 0
    end
end
