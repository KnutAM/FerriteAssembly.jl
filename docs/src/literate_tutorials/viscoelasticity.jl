using FerriteAssembly

struct ZenerMaterial{T}
    K::T    # Bulk modulus 
    G1::T   # Shear modulus, parallel
    G2::T   # Shear modulus, series 
    η::T    # Damping modulus 
end

struct ZenerState{TT<:SymmetricTensor}
    ϵv::TT # Viscous strain
end

function calculate_viscous_strain(m::ZenerMaterial, old::ZenerState, ϵ, Δt)
    return (Δt*2*m.G2*dev(ϵ) + η*old.ϵv)/(η + Δt*2*m.G2)
end

function calculate_stress(m::ZenerMaterial, old::ZenerState, ϵ, Δt)
    ϵdev = dev(ϵ)
    ϵv = calculate_viscous_strain(m, old, ϵ, Δt)
    return (m.G1+m.G2)*2*ϵdev - 2*m.G2*ϵv + m.K*vol(ϵ)
end

function FerriteAssembly.element_routine!(Ke, re, state, ae, m::ZenerMaterial, cv::CellValues, buffer)
    Δt = get_time_increment(buffer)
    old_state = get_old_state(buffer)
    for q_point in 1:getnquadpoints(cv)
        dΩ = getdetJdV(cv, q_point)
        dσdϵ, σ = gradient(e -> calculate_stress(m, old, e, Δt))
