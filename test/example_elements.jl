function get_dh_cv(grid, ::Val{false})
    dh = DofHandler(grid); add!(dh, :u, 1); close!(dh)
    cv = CellScalarValues(QuadratureRule{2,RefCube}(2), Ferrite.default_interpolation(getcelltype(grid,1)))
    return dh, cv
end
function get_dh_cv(grid, ::Val{true})
    dh = DofHandler(grid); add!(dh, :u, 2); close!(dh)
    cv = CellVectorValues(QuadratureRule{2,RefCube}(2), Ferrite.default_interpolation(getcelltype(grid,1)))
    return dh, cv
end

function assemble_test(dh, cv, m, a, aold, Δt)
    buffer, new_states, old_states = setup_assembly(dh, m, cv; a=a)
    
    # Assemble both K and r
    K = create_sparsity_pattern(dh)
    r = zeros(ndofs(dh))
    assembler = start_assemble(K, r)
    doassemble!(assembler, new_states, buffer; a=a, aold=aold, old_states=old_states, Δt=Δt)
    
    # Assemble only r
    r_direct = zeros(ndofs(dh))
    r_assembler = ReAssembler(r_direct)
    doassemble!(r_assembler, new_states, buffer; a=a, aold=aold, old_states=old_states, Δt=Δt)
    
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

# Test material from MaterialModelsBase
struct MMB_Test_Elastic{TT<:SymmetricTensor{4}} <: MMB.AbstractMaterial
    C::TT
end
MMB.material_response(m::MMB_Test_Elastic, ϵ, old, args...; kwargs...) = (m.C ⊡ ϵ, m.C, old)

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
        mmb = MMB_Test_Elastic(m.C)
        test_equality(m, weak, Val(true))
        test_equality(m, mmb, Val(true))
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
            dh = DofHandler(grid); add!(dh, :u, 2); add!(dh, :p, 1); close!(dh)
            dh_mech = DofHandler(grid); add!(dh_mech, :u, 2); close!(dh_mech)
            dh_darcy = DofHandler(grid); add!(dh_darcy, :p, 1); close!(dh_darcy)

            qr = QuadratureRule{2,RefCube}(2)
            ip = Ferrite.default_interpolation(getcelltype(grid))
            cv = (u=CellVectorValues(qr, ip), p=CellScalarValues(qr, ip))
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