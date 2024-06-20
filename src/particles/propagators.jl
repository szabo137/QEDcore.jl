import QEDbase: propagator

function _scalar_propagator(K::QEDbase.AbstractFourMomentum, mass::Real)
    return one(mass) / (K * K - mass^2)
end

function _scalar_propagator(K::QEDbase.AbstractFourMomentum)
    return one(getT(K)) / (K * K)
end

function _fermion_propagator(P::QEDbase.AbstractFourMomentum, mass::Real)
    return (slashed(P) + mass * one(DiracMatrix)) * _scalar_propagator(P, mass)
end

function _fermion_propagator(P::QEDbase.AbstractFourMomentum)
    return (slashed(P)) * _scalar_propagator(P)
end

function QEDbase.propagator(particle_type::BosonLike, K::QEDbase.AbstractFourMomentum)
    return _scalar_propagator(K, QEDbase.mass(particle_type))
end

function QEDbase.propagator(particle_type::Photon, K::QEDbase.AbstractFourMomentum)
    return _scalar_propagator(K)
end

function QEDbase.propagator(particle_type::FermionLike, P::QEDbase.AbstractFourMomentum)
    return _fermion_propagator(P, QEDbase.mass(particle_type))
end
