include("attached_values_types.jl")

# Facet routine adding the facet length to each dof attached to :u, used to check facet domains with attached variables
struct AttachedFacetTest end
function FerriteAssembly.facet_routine!(Ke, re, ae, ::AttachedFacetTest, fv, buffer)
    for q_point in 1:getnquadpoints(fv), i in dof_range(buffer, :e33)
        re[i] += getdetJdV(fv, q_point)
    end
end

# Manual assembly of generalized plane strain, where `e33` is the out-of-plane strain, for the 3d material `m`
function manual_generalized_plane_strain(dh, cv, m, a)
    K = zeros(ndofs(dh), ndofs(dh))
    r = zeros(ndofs(dh))
    ie33 = only(algebraic_dofs(dh, :e33))
    e3e3 = SymmetricTensor{2, 3}((i, j) -> (i == 3 && j == 3) ? 1.0 : 0.0)
    for cell in CellIterator(dh)
        reinit!(cv, cell)
        dofs = vcat(celldofs(cell), ie33)
        ae = a[dofs]
        n = getnbasefunctions(cv)
        B(q, i) = i ≤ n ? symmetric(lift3d(shape_gradient(cv, q, i))) : e3e3
        for q in 1:getnquadpoints(cv)
            ϵ = sum(i -> B(q, i) * ae[i], 1:(n + 1))
            σ, D, _ = MMB.material_response(m, ϵ, MMB.initial_material_state(m))
            dΩ = getdetJdV(cv, q)
            for i in 1:(n + 1)
                r[dofs[i]] += (B(q, i) ⊡ σ) * dΩ
                for j in 1:(n + 1)
                    K[dofs[i], dofs[j]] += (B(q, i) ⊡ D ⊡ B(q, j)) * dΩ
                end
            end
        end
    end
    return K, r
end

# Values type with an extra base function whose value has a different type
struct MismatchValues{CV} <: Ferrite.AbstractCellValues
    v::OutOfPlaneCellValues{CV}
end
for op in (:getnbasefunctions, :getnquadpoints, :getngeobasefunctions)
    @eval Ferrite.$op(v::MismatchValues) = Ferrite.$op(v.v)
end
Ferrite.getdetJdV(v::MismatchValues, q::Int) = getdetJdV(v.v, q)
Ferrite.geometric_value(v::MismatchValues, q::Int, i::Int) = Ferrite.geometric_value(v.v, q, i)
Ferrite.reinit!(v::MismatchValues, cell::Union{Nothing, Ferrite.AbstractCell}, x::AbstractVector{<:Vec}) = reinit!(v.v, cell, x)
Ferrite.shape_gradient(v::MismatchValues, q::Int, i::Int) = shape_gradient(v.v, q, i)
function Ferrite.shape_value(v::MismatchValues, q::Int, i::Int)
    return i < getnbasefunctions(v) ? shape_value(v.v, q, i) : zero(Vec{2})
end

# Mean value of the first component of `u` (including attached modes) enforced by the standalone variable `λ`
function add_multiplier_terms!(Ke, re, ae, cv, buffer)
    ru = dof_range(buffer, :u)
    iλ = only(dof_range(buffer, :λ))
    for q_point in 1:getnquadpoints(cv)
        dΩ = getdetJdV(cv, q_point)
        re[iλ] += function_value(cv, q_point, ae, ru)[1] * dΩ
        for (i, I) in pairs(ru)
            N = shape_value(cv, q_point, i)[1]
            re[I] += ae[iλ] * N * dΩ
            if Ke !== nothing
                Ke[I, iλ] += N * dΩ
                Ke[iλ, I] += N * dΩ
            end
        end
    end
end
struct MultiplierOnly end
FerriteAssembly.element_routine!(Ke, re, state, ae, ::MultiplierOnly, cv, buffer) = add_multiplier_terms!(Ke, re, ae, cv, buffer)
# Calls the unmodified routine of `inner` for the fields via `field_block`
struct WithMultiplier{M}
    inner::M
end
function FerriteAssembly.element_routine!(Ke, re, state, ae, m::WithMultiplier, cv, buffer)
    Ke_f, re_f, ae_f, buffer_f = FA.field_block(Ke, re, ae, buffer)
    FA.element_routine!(Ke_f, re_f, state, ae_f, m.inner, cv, buffer_f)
    add_multiplier_terms!(Ke, re, ae, cv, buffer)
end
# Only residual, calling the inner residual via `field_block`
struct WithMultiplierAD{M}
    inner::M
end
for MT in (WithMultiplier, WithMultiplierAD)
    @eval FerriteAssembly.create_cell_state(m::$MT, args...) = FerriteAssembly.create_cell_state(m.inner, args...)
    @eval FerriteAssembly.allocate_cell_cache(m::$MT, cv) = FerriteAssembly.allocate_cell_cache(m.inner, cv)
end
function FerriteAssembly.element_residual!(re, state, ae, m::WithMultiplierAD, cv, buffer)
    _, re_f, ae_f, buffer_f = FA.field_block(nothing, re, ae, buffer)
    FA.element_residual!(re_f, state, ae_f, m.inner, cv, buffer_f)
    add_multiplier_terms!(nothing, re, ae, cv, buffer)
end

# Gradient of the partner `:b`'s field, in each quadrature point
function FerriteAssembly.eval_quadpoints_cell!(vals::AbstractVector, ::Val{:partner_gradient}, state, ae, m, cv, buffer)
    pb = FA.get_coupled_buffer(buffer, :b)
    for q_point in eachindex(vals)
        vals[q_point] = function_gradient(FA.get_values(pb), q_point, FA.get_ae(pb))
    end
end

@testset "Attached algebraic variables" begin
    ip = Lagrange{RefQuadrilateral, 1}()
    qr = QuadratureRule{RefQuadrilateral}(2)
    cv_u = CellValues(qr, ip^2)
    cv_p = CellValues(qr, ip)
    xc = Vec((0.1, -0.2))
    E, ν = (2.0, 0.3)
    m_3d = EE.J2Plasticity(; E, ν, σ0 = 0.01, H = 0.2)
    m_2d = MMB.ReducedStressState(MMB.PlaneStrain(), m_3d)

    function setup_affine_dhs(; nel = (3, 2), fields = (:u => ip^2,), variables = (:H => AlgebraicVariable{Tensor{2, 2}}(), :g => AlgebraicVariable{Vec{2}}()))
        grid = generate_grid(Quadrilateral, nel)
        dh_plain = DofHandler(grid)
        dh = DofHandler(grid)
        for (name, fip) in fields
            add!(dh_plain, name, fip)
            add!(dh, name, fip)
        end
        foreach(v -> add!(dh, v.first, v.second), variables)
        close!(dh_plain)
        close!(dh)
        return dh_plain, dh
    end

    @testset "Layout" begin
        grid = generate_grid(Quadrilateral, (2, 2))
        dh = DofHandler(grid)
        add!(dh, :u, ip^2)
        add!(dh, :p, ip)
        add!(dh, :e33, AlgebraicVariable())
        add!(dh, :b, AlgebraicVariable{Vec{3}}())
        add!(dh, :A, AlgebraicVariable{SymmetricTensor{2, 3}}())
        add!(dh, :γ, AlgebraicVariable())
        close!(dh)
        sdh = only(dh.subdofhandlers)
        declarations = (
            FA.AttachAlgebraic(:A => :u; indices = ((3, 3), (2, 3), (1, 3))), :γ,
            FA.AttachAlgebraic(:e33 => :u), FA.AttachAlgebraic(:b => :p),
        )
        layout = FA.setup_layout(sdh, declarations)
        @test layout.ncelldofs == 12
        @test layout.ranges == (A = 9:11, γ = 20:20, e33 = 12:12, b = 17:19)
        # SymmetricTensor{2, 3} components: (1, 1), (2, 1), (3, 1), (2, 2), (3, 2), (3, 3)
        @test layout.components == (A = (6, 5, 3), γ = (1,), e33 = (1,), b = (1, 2, 3))
        @test layout.values.A === nothing
        @test layout.values.b isa AlgebraicValues
        @test FA.ndofs_per_item(sdh, layout) == 20
        dofrange = FA.create_dofrange(sdh, layout)
        @test dofrange == (u = 1:12, p = 13:19, A = 9:11, γ = 20:20, e33 = 12:12, b = 17:19)
        @test layout.own_ranges == (u = 1:8, p = 13:16, A = 9:11, γ = 20:20, e33 = 12:12, b = 17:19)
        @test (@inferred FA._spatial_fieldnames(keys(dofrange), layout)) == (:u, :p)
        # Attached to the last field only
        layout_p = FA.setup_layout(sdh, (FA.AttachAlgebraic(:b => :p; indices = (3, 1)),))
        @test FA.create_dofrange(sdh, layout_p) == (u = 1:8, p = 9:14, b = 13:14)
        @test layout_p.components == (b = (3, 1),)
        # Without attached variables
        layout_γ = FA.setup_layout(sdh, (:γ,))
        @test layout_γ.fields === nothing
        @test FA.create_dofrange(sdh, layout_γ) == (u = 1:8, p = 9:12, γ = 13:13)

        # Filling the local dofs
        d = FA.DofBuffer(sdh, layout)
        @test length(d.scratch) == 12
        @test isempty(FA.DofBuffer(sdh, layout_γ).scratch)
        @test isempty(FA.DofBuffer(sdh).scratch)
        for cellnr in 1:getncells(grid)
            FA.set_dofs!(d, dh, cellnr)
            cdofs = celldofs(dh, cellnr)
            @test d.dofs[1:8] == cdofs[dof_range(sdh, :u)]
            @test d.dofs[13:16] == cdofs[dof_range(sdh, :p)]
            @test d.dofs[9:11] == algebraic_dofs(dh, :A)[[6, 5, 3]]
            @test d.dofs[12:12] == algebraic_dofs(dh, :e33)
            @test d.dofs[17:19] == algebraic_dofs(dh, :b)
            @test d.dofs[20:20] == algebraic_dofs(dh, :γ)
        end
        @test (@allocated FA.set_dofs!(d, dh, 1)) == 0
        d_p = FA.DofBuffer(sdh, layout_p)
        FA.set_dofs!(d_p, dh, 3)
        @test d_p.dofs == vcat(celldofs(dh, 3), algebraic_dofs(dh, :b)[[3, 1]])
    end

    @testset "Item buffers" begin
        dh_plain, dh = setup_affine_dhs(; fields = (:u => ip^2, :p => ip), variables = (:H => AlgebraicVariable{Tensor{2, 2}}(), :γ => AlgebraicVariable()))
        H_indices = ((2, 1), (1, 1))
        values = (u = AffineCellValues(cv_u, xc, H_indices), p = cv_p)
        declarations = (:γ, FA.AttachAlgebraic(:H => :u; indices = H_indices))
        a = rand(ndofs(dh))
        aold = rand(ndofs(dh))
        for threading in (false, true)
            db = setup_domainbuffer(DomainSpec(dh, nothing, values; algebraic_variables = declarations); threading)
            buffer = FA.get_base(FA.get_itembuffer(db))
            @test Ferrite.getfieldnames(buffer) == (:u, :p, :γ, :H)
            @test FA._spatial_fieldnames(buffer) == (:u, :p)
            @test dof_range(buffer, :u) == 1:10
            @test dof_range(buffer, :H) == 9:10
            @test dof_range(buffer, :p) == 11:14
            @test dof_range(buffer, :γ) == 15:15
            @test FA.field_dof_range(buffer, :u) == 1:8
            @test FA.field_dof_range(buffer, :p) == 11:14
            @test FA.field_dof_range(buffer, :H) == 9:10
            @test_throws ArgumentError FA.get_algebraic_values(buffer, :H)
            @test FA.get_algebraic_values(buffer, :γ) isa AlgebraicValues
            sim = Simulation(db, a, aold)
            FA.refresh_algebraic_dofs!(sim)
            for cellnr in (1, 4)
                FA.reinit_buffer!(buffer, sim, cellnr)
                cdofs = celldofs(dh, cellnr)
                H_dofs = algebraic_dofs(dh, :H)[[2, 1]] # Tensor{2,2} components: (1, 1), (2, 1), (1, 2), (2, 2)
                expected = vcat(cdofs[1:8], H_dofs, cdofs[9:12], algebraic_dofs(dh, :γ))
                @test celldofs(buffer) == expected
                @test FA.get_ae(buffer) == a[expected]
                @test FA.get_aeold(buffer) == aold[expected]
            end
            alloc_field_dof_range(b) = @allocated FA.field_dof_range(b, :u)
            alloc_field_dof_range(buffer)
            @test alloc_field_dof_range(buffer) == 0
        end
        # Without layout
        buffer = FA.get_itembuffer(setup_domainbuffer(DomainSpec(dh_plain, nothing, (u = cv_u, p = cv_p))))
        @test FA.field_dof_range(buffer, :p) == dof_range(buffer, :p) == 9:12
    end

    @testset "Validation" begin
        dh_plain, dh = setup_affine_dhs(; fields = (:u => ip^2, :p => ip), variables = (:H => AlgebraicVariable{Tensor{2, 2}}(), :γ => AlgebraicVariable(), :S => AlgebraicVariable{SymmetricTensor{2, 2}}()))
        H_indices = ((2, 1), (1, 1))
        values = (u = AffineCellValues(cv_u, xc, H_indices), p = cv_p)
        attach_H(; indices = H_indices) = FA.AttachAlgebraic(:H => :u; indices)
        spec(decl; vals = values) = DomainSpec(dh, nothing, vals; algebraic_variables = decl)
        @test setup_domainbuffer(spec((attach_H(),))) isa FA.DomainBuffer
        # Declarations
        @test_throws ArgumentError FA.AttachAlgebraic(:H => :u; indices = [1, 2])
        @test_throws ArgumentError spec((:H => :u,))
        @test_throws ArgumentError setup_domainbuffer(spec((FA.AttachAlgebraic(:H => :v; indices = H_indices),)))
        @test_throws ArgumentError setup_domainbuffer(spec((FA.AttachAlgebraic(:u => :p),)))
        @test_throws ArgumentError setup_domainbuffer(spec((FA.AttachAlgebraic(:Q => :u),)))
        @test_throws ArgumentError setup_domainbuffer(spec((attach_H(; indices = ()),)))
        @test_throws ArgumentError setup_domainbuffer(spec((attach_H(; indices = ((3, 1), (1, 1))),)))
        @test_throws ArgumentError setup_domainbuffer(spec((attach_H(; indices = ((1, 1), (1, 1))),)))
        @test_throws ArgumentError setup_domainbuffer(spec((attach_H(; indices = (1, 2)),)))
        @test_throws ArgumentError setup_domainbuffer(spec((attach_H(; indices = ((1, 1, 1), (2, 1))),)))
        @test_throws ArgumentError setup_domainbuffer(spec((FA.AttachAlgebraic(:S => :u; indices = ((1, 2), (2, 1))),)))
        @test_throws ArgumentError setup_domainbuffer(spec((FA.AttachAlgebraic(:γ => :u; indices = (2,)),)))
        # A variable once per domain (also across modes and fields)
        @test_throws ArgumentError setup_domainbuffer(spec((attach_H(), :H)))
        @test_throws ArgumentError setup_domainbuffer(spec((attach_H(), FA.AttachAlgebraic(:H => :p))))
        # Values
        @test_throws ArgumentError setup_domainbuffer(spec((attach_H(),); vals = (u = cv_u, p = cv_p)))
        @test_throws ArgumentError setup_domainbuffer(spec((attach_H(),); vals = (v = values.u, p = cv_p)))
        @test_throws ArgumentError setup_domainbuffer(spec((attach_H(),); vals = values.u))
        # Single field values need not be NamedTuples
        db = setup_domainbuffer(DomainSpec(setup_affine_dhs()[2], nothing, values.u; algebraic_variables = (attach_H(),)))
        @test dof_range(FA.get_itembuffer(db), :u) == 1:10
        # Values of fields without attached variables are not checked
        @test setup_domainbuffer(spec((:γ,); vals = (x = cv_u, y = cv_p))) isa FA.DomainBuffer
    end

    @testset "Generalized plane strain, unmodified routine" begin
        grid = generate_grid(Quadrilateral, (3, 2), Vec((0.0, 0.0)), Vec((2.0, 1.0)))
        dh = DofHandler(grid)
        add!(dh, :u, ip^2)
        add!(dh, :e33, AlgebraicVariable())
        close!(dh)
        m_elastic = EE.J2Plasticity(; E, ν, σ0 = Inf, H = 1.0)
        cv_oop = OutOfPlaneCellValues(CellValues(qr, ip^2))
        for m in (m_elastic, m_3d)
            db = setup_domainbuffer(DomainSpec(dh, m, cv_oop; algebraic_variables = (FA.AttachAlgebraic(:e33 => :u),)))
            K = allocate_matrix(dh; algebraic_couplings = FA.algebraic_couplings(db))
            r = zeros(ndofs(dh))
            a = 1.0e-3 * randn(ndofs(dh))
            work!(KeReAssembler(K, r), db; a)
            K_ref, r_ref = manual_generalized_plane_strain(dh, CellValues(qr, ip^2), m, a)
            @test K ≈ K_ref
            @test r ≈ r_ref
        end

        # Homogeneous strain: the out-of-plane strain gives zero out-of-plane stress
        db = setup_domainbuffer(DomainSpec(dh, m_elastic, cv_oop; algebraic_variables = (FA.AttachAlgebraic(:e33 => :u),)))
        ch = ConstraintHandler(dh)
        ϵ1, ϵ2 = (0.01, -0.004)
        add!(ch, Dirichlet(:u, union(getfacetset.((grid,), ("left", "right", "top", "bottom"))...), x -> Vec((ϵ1 * x[1], ϵ2 * x[2]))))
        close!(ch)
        update!(ch, 0.0)
        K = allocate_matrix(dh, ch; algebraic_couplings = FA.algebraic_couplings(db))
        r = zeros(ndofs(dh))
        a = zeros(ndofs(dh))
        apply!(a, ch)
        work!(KeReAssembler(K, r), db; a)
        apply_zero!(K, r, ch)
        a .-= K \ r
        @test a[only(algebraic_dofs(dh, :e33))] ≈ -ν / (1 - ν) * (ϵ1 + ϵ2)

        # Residual scaling uses the field's own dofs
        scaling = ElementResidualScaling(dh, 1)
        work!(KeReAssembler(K, r; scaling), db; a)
        factors = Dict(:u => 0.0, :e33 => 0.0)
        buffer = FA.get_itembuffer(db)
        sim = Simulation(db, a)
        for cellnr in 1:getncells(grid)
            FA.reinit_buffer!(buffer, sim, cellnr)
            FA.element_routine!(FA.get_Ke(buffer), FA.get_re(buffer), FA.get_state(buffer), FA.get_ae(buffer), m_elastic, cv_oop, buffer)
            factors[:u] += sum(abs, FA.get_re(buffer)[1:8])
            factors[:e33] += abs(FA.get_re(buffer)[9])
        end
        @test scaling.factors[:u] ≈ factors[:u]
        @test scaling.factors[:e33] ≈ factors[:e33]
    end

    @testset "Affine modes, unmodified routines" begin
        H_indices = ((2, 1), (1, 1), (1, 2)) # Reordered subset
        g_indices = (2,)
        dh_plain, dh = setup_affine_dhs()
        T = affine_map(dh_plain, dh, xc, H_indices, g_indices)
        declarations = (FA.AttachAlgebraic(:H => :u; indices = H_indices), FA.AttachAlgebraic(:g => :u; indices = g_indices))
        values = AffineCellValues(cv_u, xc, H_indices, g_indices)
        weakform = EE.WeakForm((δu, ∇δu, u, ∇u, u_dot, ∇u_dot) -> (∇δu ⊡ (2 * symmetric(∇u) + tr(∇u) * one(∇u)) + δu ⋅ u_dot))
        a = 1.0e-3 * randn(ndofs(dh))
        aold = 1.0e-3 * randn(ndofs(dh))
        for (material, autodiffbuffer, threading) in ((m_2d, false, false), (m_2d, false, true), (weakform, false, false), (weakform, true, false), (weakform, true, true))
            db_plain = setup_domainbuffer(DomainSpec(dh_plain, material, cv_u); autodiffbuffer)
            Kp = allocate_matrix(dh_plain)
            rp = zeros(ndofs(dh_plain))
            set_time_increment!(db_plain, 0.1)
            work!(KeReAssembler(Kp, rp), db_plain; a = T * a, aold = T * aold)

            db = setup_domainbuffer(DomainSpec(dh, material, values; algebraic_variables = declarations); autodiffbuffer, threading)
            K = allocate_matrix(dh; algebraic_couplings = FA.algebraic_couplings(db))
            r = zeros(ndofs(dh))
            set_time_increment!(db, 0.1)
            work!(KeReAssembler(K, r; atomic = threading), db; a, aold)
            @test K ≈ T' * Kp * T
            @test r ≈ T' * rp
        end

        # Multiple fields, attached to the first field (the routine for PoroElasticPlaneStrain uses AD)
        dh_plain, dh = setup_affine_dhs(; fields = (:u => ip^2, :p => ip))
        T = affine_map(dh_plain, dh, xc, H_indices, g_indices)
        poro = EE.PoroElasticPlaneStrain(; E, ν, k = 0.5, α = 0.7, β = 0.3)
        a = 1.0e-3 * randn(ndofs(dh))
        aold = 1.0e-3 * randn(ndofs(dh))
        for autodiffbuffer in (false, true)
            db_plain = setup_domainbuffer(DomainSpec(dh_plain, poro, (u = cv_u, p = cv_p)); autodiffbuffer)
            Kp = allocate_matrix(dh_plain)
            rp = zeros(ndofs(dh_plain))
            set_time_increment!(db_plain, 0.1)
            work!(KeReAssembler(Kp, rp), db_plain; a = T * a, aold = T * aold)

            db = setup_domainbuffer(DomainSpec(dh, poro, (u = values, p = cv_p); algebraic_variables = declarations); autodiffbuffer)
            K = allocate_matrix(dh; algebraic_couplings = FA.algebraic_couplings(db))
            r = zeros(ndofs(dh))
            set_time_increment!(db, 0.1)
            work!(KeReAssembler(K, r), db; a, aold)
            @test K ≈ T' * Kp * T
            @test r ≈ T' * rp
        end
    end

    @testset "Attached and standalone variables (field_block)" begin
        H_indices = ((2, 1), (1, 1))
        dh_plain, dh = setup_affine_dhs(; variables = (:H => AlgebraicVariable{Tensor{2, 2}}(), :λ => AlgebraicVariable()))
        values = AffineCellValues(cv_u, xc, H_indices)
        attach_H = FA.AttachAlgebraic(:H => :u; indices = H_indices)
        weakform = EE.WeakForm((δu, ∇δu, u, ∇u, u_dot, ∇u_dot) -> (∇δu ⊡ (2 * symmetric(∇u) + tr(∇u) * one(∇u)) + δu ⋅ u_dot))
        a = 1.0e-3 * randn(ndofs(dh))
        aold = 1.0e-3 * randn(ndofs(dh))
        function assemble(material, declarations; kwargs...)
            db = setup_domainbuffer(DomainSpec(dh, material, values; algebraic_variables = declarations); kwargs...)
            K = allocate_matrix(dh; algebraic_couplings = FA.algebraic_couplings(db))
            r = zeros(ndofs(dh))
            set_time_increment!(db, 0.1)
            work!(KeReAssembler(K, r; atomic = true), db; a, aold)
            return K, r
        end
        K_mult, r_mult = assemble(MultiplierOnly(), (attach_H, :λ))
        for inner in (m_2d, weakform)
            K_inner, r_inner = assemble(inner, (attach_H,))
            for (outer, autodiffbuffer, threading) in (
                    (WithMultiplier(inner), false, false), (WithMultiplier(inner), true, false),
                    (WithMultiplier(inner), false, true), (WithMultiplierAD(inner), true, false),
                    (WithMultiplierAD(inner), false, true),
                )
                K, r = assemble(outer, (attach_H, :λ); autodiffbuffer, threading)
                @test K ≈ K_inner + K_mult
                @test r ≈ r_inner + r_mult
            end
        end

        # The buffer view
        db = setup_domainbuffer(DomainSpec(dh, nothing, values; algebraic_variables = (:λ, attach_H)))
        buffer = FA.get_itembuffer(db)
        sim = Simulation(db, a, aold)
        FA.refresh_algebraic_dofs!(sim)
        FA.reinit_buffer!(buffer, sim, 2)
        Ke, re, ae = (FA.get_Ke(buffer), FA.get_re(buffer), FA.get_ae(buffer))
        Ke_f, re_f, ae_f, buffer_f = FA.field_block(Ke, re, ae, buffer)
        @test size(Ke_f) == (10, 10)
        @test length(re_f) == length(ae_f) == 10
        @test ae_f == ae[1:10]
        @test FA.get_ae(buffer_f) == FA.get_ae(buffer)[1:10]
        @test FA.get_aeold(buffer_f) == FA.get_aeold(buffer)[1:10]
        @test celldofs(buffer_f) == celldofs(buffer)[1:10]
        @test size(FA.get_Ke(buffer_f)) == (10, 10)
        @test length(FA.get_re(buffer_f)) == 10
        @test dof_range(buffer_f, :u) == 1:10
        @test FA.get_values(buffer_f) === FA.get_values(buffer)
        @test FA.get_material(buffer_f) === FA.get_material(buffer)
        @test Ferrite.cellid(buffer_f) == 2
        @test getcoordinates(buffer_f) === getcoordinates(buffer)
        alloc_field_block(Ke, re, ae, buffer) = @allocated FA.field_block(Ke, re, ae, buffer)
        alloc_field_block(Ke, re, ae, buffer)
        @test alloc_field_block(Ke, re, ae, buffer) == 0
        # Without standalone variables, and without algebraic variables
        for decl in ((attach_H,), ())
            db = setup_domainbuffer(DomainSpec(decl == () ? dh_plain : dh, nothing, decl == () ? cv_u : values; algebraic_variables = decl))
            buffer = FA.get_itembuffer(db)
            n = length(FA.get_ae(buffer))
            Ke_f, re_f, ae_f, buffer_f = FA.field_block(FA.get_Ke(buffer), FA.get_re(buffer), FA.get_ae(buffer), buffer)
            @test size(Ke_f) == (n, n)
            @test length(FA.get_ae(buffer_f)) == n
        end
    end

    @testset "Integrators" begin
        H_indices = ((1, 2),)
        dh_plain, dh = setup_affine_dhs(; variables = (:H => AlgebraicVariable{Tensor{2, 2}}(),))
        T = affine_map(dh_plain, dh, xc, H_indices)
        a = randn(ndofs(dh))
        values = AffineCellValues(cv_u, xc, H_indices)
        db = setup_domainbuffer(DomainSpec(dh, nothing, values; algebraic_variables = (FA.AttachAlgebraic(:H => :u; indices = H_indices),)))
        db_plain = setup_domainbuffer(DomainSpec(dh_plain, nothing, cv_u))
        fun = (u, ∇u, state) -> (u, ∇u)
        val0 = (zero(Vec{2}), zero(Tensor{2, 2}))
        integrator = SimpleIntegrator(fun, val0)
        work!(integrator, db; a)
        integrator_plain = SimpleIntegrator(fun, val0)
        work!(integrator_plain, db_plain; a = T * a)
        @test integrator.val[1] ≈ integrator_plain.val[1]
        @test integrator.val[2] ≈ integrator_plain.val[2]
        qe = QuadPointEvaluator{Float64}(db, (m, u, ∇u, state) -> u[1] + ∇u[1, 2])
        work!(qe, db; a)
        qe_plain = QuadPointEvaluator{Float64}(db_plain, (m, u, ∇u, state) -> u[1] + ∇u[1, 2])
        work!(qe_plain, db_plain; a = T * a)
        @test qe.data ≈ qe_plain.data
    end

    @testset "Coupled simulations" begin
        H_indices = ((1, 2), (2, 2))
        dh_plain, dh = setup_affine_dhs(; variables = (:H => AlgebraicVariable{Tensor{2, 2}}(),))
        T = affine_map(dh_plain, dh, xc, H_indices)
        dh_reader = DofHandler(dh.grid)
        add!(dh_reader, :c, ip)
        close!(dh_reader)
        a = randn(ndofs(dh))
        for threading in (false, true), autodiffbuffer in (false, true)
            partner = DomainSpec(dh, m_2d, AffineCellValues(cv_u, xc, H_indices); algebraic_variables = (FA.AttachAlgebraic(:H => :u; indices = H_indices),))
            results = map(((partner, a), (DomainSpec(dh_plain, m_2d, cv_u), T * a))) do (spec, a_partner)
                sim_reader = Simulation(setup_domainbuffer(DomainSpec(dh_reader, nothing, cv_p); threading, autodiffbuffer))
                sim_partner = Simulation(setup_domainbuffer(spec; threading, autodiffbuffer), a_partner)
                group = CoupledSimulations((r = sim_reader,); refs = (b = sim_partner,))
                qe = QuadPointEvaluator{Tensor{2, 2, Float64, 4}}(group.r.db, :partner_gradient)
                work!(qe, group.r)
                return qe.data
            end
            @test all(cellnr -> results[1][cellnr] ≈ results[2][cellnr], 1:getncells(dh.grid))
        end
    end

    @testset "Renumbering and replacing material" begin
        H_indices = ((2, 2), (1, 2))
        dh_plain, dh = setup_affine_dhs(; variables = (:H => AlgebraicVariable{Tensor{2, 2}}(),))
        declarations = (FA.AttachAlgebraic(:H => :u; indices = H_indices),)
        values = AffineCellValues(cv_u, xc, H_indices)
        db = setup_domainbuffer(DomainSpec(dh, m_2d, values; algebraic_variables = declarations))
        a = av_renumber!(dh, 1.0e-3 * randn(ndofs(dh)))
        db_new = setup_domainbuffer(DomainSpec(dh, m_2d, values; algebraic_variables = declarations))
        K1, K2 = (allocate_matrix(dh; algebraic_couplings = FA.algebraic_couplings(db_new)) for _ in 1:2)
        r1, r2 = (zeros(ndofs(dh)) for _ in 1:2)
        work!(KeReAssembler(K1, r1), db; a)
        work!(KeReAssembler(K2, r2), db_new; a)
        @test K1 ≈ K2
        @test r1 ≈ r2
        m_elastic = MMB.ReducedStressState(MMB.PlaneStrain(), EE.J2Plasticity(; E, ν, σ0 = Inf, H = 1.0))
        db_elastic = FA.replace_material(db, _ -> m_elastic)
        @test dof_range(FA.get_itembuffer(db_elastic), :u) == 1:10
        work!(KeReAssembler(K1, r1), db_elastic; a)
        @test K1 ≈ K2
    end

    @testset "Facet domain" begin
        grid = generate_grid(Quadrilateral, (3, 2), Vec((0.0, 0.0)), Vec((2.0, 1.0)))
        dh = DofHandler(grid)
        add!(dh, :u, ip^2)
        add!(dh, :e33, AlgebraicVariable())
        close!(dh)
        fv = FacetValues(FacetQuadratureRule{RefQuadrilateral}(2), ip^2)
        spec(vals) = DomainSpec(dh, AttachedFacetTest(), vals; set = getfacetset(grid, "right"), algebraic_variables = (FA.AttachAlgebraic(:e33 => :u),))
        @test_throws ArgumentError setup_domainbuffer(spec(fv))
        db = setup_domainbuffer(spec(PaddedFacetValues(fv, 1)))
        buffer = FA.get_itembuffer(db)
        @test dof_range(buffer, :u) == 1:9
        @test dof_range(buffer, :e33) == 9:9
        K = allocate_matrix(dh; algebraic_couplings = FA.algebraic_couplings(db))
        r = zeros(ndofs(dh))
        work!(KeReAssembler(K, r), db)
        @test r[only(algebraic_dofs(dh, :e33))] ≈ 1.0
        @test celldofs(buffer) == vcat(celldofs(dh, FA.cellid(buffer)), algebraic_dofs(dh, :e33))
    end

    @testset "Values interface" begin
        grid = generate_grid(Quadrilateral, (1, 1))
        x = getcoordinates(grid, 1)
        cell = getcells(grid, 1)
        @test FA.check_values_interface(CellValues(qr, ip), cell, x)
        for base_ip in (ip, ip^2)
            cvh = CellValues(qr, base_ip; update_hessians = true)
            @test FA.check_values_interface(OutOfPlaneCellValues(cvh), cell, x; hessians = true)
        end
        cvh = CellValues(qr, ip^2; update_hessians = true)
        @test FA.check_values_interface(AffineCellValues(cvh, xc, ((2, 1), (1, 1)), (2,)), cell, x; hessians = true)
        @test_throws ErrorException FA.check_values_interface(MismatchValues(OutOfPlaneCellValues(CellValues(qr, ip^2))), cell, x)
    end
end
