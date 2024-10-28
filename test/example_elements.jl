function get_dh_cv(grid, ::Val{false})
    ip = geometric_interpolation(getcelltype(grid,1))
    dh = DofHandler(grid); add!(dh, :u, ip); close!(dh)
    RefShape = Ferrite.getrefshape(ip)
    cv = CellValues(QuadratureRule{RefShape}(2), ip)
    return dh, cv
end
function get_dh_cv(grid, ::Val{true})
    ip = VectorizedInterpolation(geometric_interpolation(getcelltype(grid,1)))
    dh = DofHandler(grid); add!(dh, :u, ip); close!(dh)
    RefShape = Ferrite.getrefshape(ip)
    cv = CellValues(QuadratureRule{RefShape}(2), ip)
    return dh, cv
end

function assemble_test(dh, cv, m, a, aold, Δt)
    buffer = setup_domainbuffer(DomainSpec(dh, m, cv); a=a)
    set_time_increment!(buffer, Δt)
    # Assemble both K and r
    K = allocate_matrix(dh)
    r = zeros(ndofs(dh))
    assembler = start_assemble(K, r)
    work!(assembler, buffer; a=a, aold=aold)
    
    # Assemble only r
    r_direct = zeros(ndofs(dh))
    r_assembler = ReAssembler(r_direct)
    work!(r_assembler, buffer; a=a, aold=aold)
    
    return K, r, r_direct
end

function test_equality(m1, m2, is_vector::Val)
    grid = generate_grid(Quadrilateral, (2,3))
    dh, cv = get_dh_cv(grid, is_vector)
    a = rand(ndofs(dh))
    aold = rand(ndofs(dh))
    Δt = rand()
    K1, r1, r1_direct = assemble_test(dh, cv, m1, a, aold, Δt)
    K2, r2, r2_direct = assemble_test(dh, cv, m2, a, aold, Δt)
    
    # Check equivalence between assembling stiffness and residual vs just residual
    @test r1 ≈ r1_direct # (element_routine! vs element_residual!)
    @test r2 ≈ r2_direct # (in some cases, the same routine is actually called.)
    @test r1 ≈ r2
    @test K1 ≈ K2
end

@testset "HeatEquation" begin
    @testset "StationaryFourier" begin
        k = rand()
        m = EE.StationaryFourier(k)
        weak = EE.WeakForm((δu, ∇δu, u, ∇u, u_dot, ∇u_dot) -> k*(∇δu ⋅ ∇u))
        test_equality(m, weak, Val(false))
    end
    @testset "TransientFourier" begin
        k, c = rand(2)
        m = EE.TransientFourier(k, c)
        weak = EE.WeakForm((δu, ∇δu, u, ∇u, u_dot, ∇u_dot) -> δu*c*u_dot + k*(∇δu ⋅ ∇u))
        test_equality(m, weak, Val(false))
    end
end

@testset "Mechanical" begin
    @testset "LinearElasticity" begin
        E = 1.0 + rand()
        ν = 0.1 + rand()/3 #∈ [0.1, 0.433]
        G = E/(2*(1+ν)); K = E/(3*(1-2ν))
        m = EE.ElasticPlaneStrain(;E=E, ν=ν)
        weak = EE.WeakForm((δu, ∇δu, u, ∇u, u_dot, ∇u_dot) -> (∇δu ⊡ (2*G*dev(symmetric(∇u)) + 3*K*vol(∇u))))
        mmb = MMB.ReducedStressState(MMB.PlaneStrain(),
            EE.J2Plasticity(;E=E, ν=ν, σ0=Inf, H=1.0))
        test_equality(m, weak, Val(true))
        test_equality(m, mmb, Val(true))
    end
    @testset "J2Plasticity" begin
        E = 1.0 + rand()
        H = 0.1 + 0.9*rand()
        ν = 0.1 + rand()/3 #∈ [0.1, 0.433]
        σ0 = E/(100*(1+rand()))
        m = MMB.ReducedStressState(MMB.UniaxialStress(),
            EE.J2Plasticity(;E, ν, σ0, H))

        ϵv = collect(range(0, 2σ0/E, 100))
        σv = Float64[]
        dσdϵ_v = Float64[]
        ϵ22v = Float64[]
        state = MMB.initial_material_state(m)
        for ϵ in ϵv 
            ϵt = SymmetricTensor{2,1}((ϵ,))
            σ, dσdϵ, state, ϵf = MMB.material_response(m, ϵt, state)
            push!(σv, σ[1,1])
            push!(dσdϵ_v, dσdϵ[1,1,1,1])
            push!(ϵ22v, ϵf[2,2])
        end
        Hp = E*H/(E+H)
        @test dσdϵ_v[1] ≈ E 
        @test σv[2]-σv[1] ≈ E*(ϵv[2]-ϵv[1])
        @test ϵ22v[2] ≈ -ν*ϵv[2]
        @test dσdϵ_v[end] ≈ Hp
        @test σv[end]-σv[end-1] ≈ Hp*(ϵv[end]-ϵv[end-1])
    end
end

@testset "PorousMedia" begin
    @testset "PoroElasticPlaneStrain" begin
        E = rand() + 1.0
        ν = rand()/3 + 0.1
        k, β = rand(2)
        for α in (0.0, rand()) # α = 0.0 => No coupling => correct residual values
            m = EE.PoroElasticPlaneStrain(;E=E, ν=ν, k=k, α=α, β=β)
            m_mech = EE.ElasticPlaneStrain(;E=E, ν=ν)
            m_darcy = EE.WeakForm((δp, ∇δp, p, ∇p, p_dot, ∇p_dot) -> δp*β*p_dot + k*(∇δp ⋅ ∇p))

            grid = generate_grid(Quadrilateral, (1,1)) # Single element test
            ip = Lagrange{RefQuadrilateral,1}()
            dh = DofHandler(grid); add!(dh, :u, ip^2); add!(dh, :p, ip); close!(dh)
            dh_mech = DofHandler(grid); add!(dh_mech, :u, ip^2); close!(dh_mech)
            dh_darcy = DofHandler(grid); add!(dh_darcy, :p, ip); close!(dh_darcy)

            qr = QuadratureRule{Ferrite.getrefshape(ip)}(2)
            cv = (u=CellValues(qr, ip^2), p=CellValues(qr, ip))
            i_mech = 1:ndofs(dh_mech); i_darcy = ndofs(dh_mech) .+ (1:ndofs(dh_darcy))
            a = rand(ndofs(dh)); a_mech = a[i_mech]; a_darcy = a[i_darcy]
            aold = rand(ndofs(dh)); aold_mech = aold[i_mech]; aold_darcy = aold[i_darcy];
            Δt = rand()
            K, r = assemble_test(dh, cv, m, a, aold, Δt)
            K_mech, r_mech = assemble_test(dh_mech, cv[:u], m_mech, a_mech, aold_mech, Δt)
            K_darcy, r_darcy = assemble_test(dh_darcy, cv[:p], m_darcy, a_darcy, aold_darcy, Δt)
            @test K[i_mech, i_mech] ≈ K_mech 
            @test K[i_darcy, i_darcy] ≈ K_darcy 
            if α == 0.0 # No coupling
                # Should be no coupling, so check this first
                @test norm(K[i_mech,i_darcy]) < 1e-10*norm(K)
                @test norm(K[i_darcy,i_mech]) < 1e-10*norm(K)
                # As there is no coupling, the residuals should be equal
                @test r[i_mech] ≈ r_mech
                @test r[i_darcy] ≈ r_darcy                
            end
        end
    end
end