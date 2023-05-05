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
    K = create_sparsity_pattern(dh)
    r = zeros(ndofs(dh))
    states = create_states(dh, m, cv, a)
    cellbuffer = setup_cellbuffer(dh, cv, m)
    assembler = start_assemble(K,r)
    doassemble!(assembler, cellbuffer, states, dh, a, aold, Δt);
    return K, r
end

function test_equality(m1, m2, is_vector::Val)
    grid = generate_grid(Quadrilateral, (2,3))
    dh, cv = get_dh_cv(grid, is_vector)
    a = rand(ndofs(dh))
    aold = rand(ndofs(dh))
    Δt = rand()
    K1, r1 = assemble_test(dh, cv, m1, a, aold, Δt)
    K2, r2 = assemble_test(dh, cv, m2, a, aold, Δt)
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
        weak = EE.WeakForm((δu, ∇δu, u, ∇u, u_dot, ∇u_dot) -> -k*(∇δu ⋅ ∇u))
        test_equality(m, weak, Val(false))
    end
    @testset "TransientFourier" begin
        k, c = rand(2)
        m = EE.TransientFourier(k, c)
        weak = EE.WeakForm((δu, ∇δu, u, ∇u, u_dot, ∇u_dot) -> δu*c*u_dot - k*(∇δu ⋅ ∇u))
        test_equality(m, weak, Val(false))
    end
end

@testset "Mechanical" begin
    @testset "LinearElasticity" begin
        E = 1.0 + rand()
        ν = 0.1 + rand()/3 #∈ [0.1, 0.433]
        G = E/(2*(1+ν)); K = E/(3*(1-2ν))
        m = EE.ElasticPlaneStress(;E=E, ν=ν)
        weak = EE.WeakForm((δu, ∇δu, u, ∇u, u_dot, ∇u_dot) -> (∇δu ⊡ (2*G*dev(symmetric(∇u)) + 3*K*vol(∇u))))
        mmb = MMB_Test_Elastic(m.C)
        test_equality(m, weak, Val(true))
        test_equality(m, mmb, Val(true))
    end
end