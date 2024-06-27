function _scalar_propagator(K::AbstractFourMomentum, mass::Real)
    return one(mass) / (K * K - mass^2)
end

function _scalar_propagator(K::AbstractFourMomentum)
    return one(getT(K)) / (K * K)
end

function _fermion_propagator(P::AbstractFourMomentum, mass::Real)
    return (slashed(P) + mass * one(DiracMatrix)) * _scalar_propagator(P, mass)
end

function _fermion_propagator(P::AbstractFourMomentum)
    return (slashed(P)) * _scalar_propagator(P)
end

function QEDbase.propagator(particle_type::BosonLike, K::AbstractFourMomentum)
    return _scalar_propagator(K, mass(particle_type))
end

function QEDbase.propagator(particle_type::Photon, K::AbstractFourMomentum)
    return _scalar_propagator(K)
end

function QEDbase.propagator(particle_type::FermionLike, P::AbstractFourMomentum)
    return _fermion_propagator(P, mass(particle_type))
end
