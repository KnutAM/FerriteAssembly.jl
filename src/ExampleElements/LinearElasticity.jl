struct ElasticPlaneStress{T} <: AbstractMaterial
    dσdϵ::SymmetricTensor{4,2,T,9}
end
function ElasticPlaneStress(;E=2.e3, ν=0.3)
    G = E / 2(1 + ν)
    K = E / 3(1 - 2ν)
    I2 = one(SymmetricTensor{2,2})
    I4vol = I2⊗I2
    I4dev = minorsymmetric(otimesu(I2,I2)) - I4vol / 3
    return ElasticPlaneStress(2G*I4dev + K*I4vol)
end;

function FerriteAssembly.element_routine!(Ke, re, state, ae, material::ElasticPlaneStress, cv::CellVectorValues, dh_fh, Δt, buffer)
    for q_point in 1:getnquadpoints(cv_u)
        dΩ = getdetJdV(cv, q_point)
        ϵ = function_symmetric_gradient(cv_u, q_point, ae)
        σ = material.dσdϵ ⊡ ϵ
        ## Assemble residual contributions
        for i in 1:getnbasefunctions(cv)
            ∇δNu = shape_symmetric_gradient(cv_u, q_point, i)
            re[i] += (∇δNu ⊡ σ )*dΩ
            ∇δNu_dσdϵ = ∇δNu ⊡  material.dσdϵ
            for j in 1:getnbasefunctions(cv)
                ∇Nu = shape_symmetric_gradient(cv_u, q_point, j)
                Ke[j,i] += (∇δNu_dσdϵ ⊡ ∇Nu)*dΩ # Since Ke is symmetric, we calculate Ke' to index faster
            end
        end
    end
end