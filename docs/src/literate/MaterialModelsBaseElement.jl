# ## Element routine
# The element routine for any material that follows the `MaterialModelsBase.jl` interface 
# can be coded as follows:
function FerriteAssembly.element_routine!(
    Ke::AbstractMatrix, re::AbstractVector, state::Vector{<:AbstractMaterialState}, 
    ae::AbstractVector, material::AbstractMaterial, cellvalues::CellVectorValues, dh_fh, Δt, buffer::CellBuffer)
    cache = buffer.cache
    n_basefuncs = getnbasefunctions(cellvalues)
    for q_point in 1:getnquadpoints(cellvalues)
        ## For each integration point, compute stress and material stiffness
        ϵ = function_symmetric_gradient(cellvalues, q_point, ae) # Total strain
        σ, D, state[q_point] = material_response(material, ϵ, state[q_point], Δt, cache)

        dΩ = getdetJdV(cellvalues, q_point)
        for i in 1:n_basefuncs
            δϵ = shape_symmetric_gradient(cellvalues, q_point, i)
            re[i] += (δϵ ⊡ σ) * dΩ ## add internal force to residual
            for j in 1:n_basefuncs
                Δϵ = shape_symmetric_gradient(cellvalues, q_point, j)
                Ke[i, j] += δϵ ⊡ D ⊡ Δϵ * dΩ
            end
        end
    end
end;
# Note that the boundary traction was not included, as this can be handled separately
# with [FerriteNeumann.jl](https://github.com/KnutAM/FerriteNeumann.jl)
