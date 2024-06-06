function _booster_fermion(mom::QEDbase.AbstractFourMomentum, mass::Real)
    return (slashed(mom) + mass * one(DiracMatrix)) / (sqrt(abs(getT(mom)) + mass))
end

function _booster_antifermion(mom::QEDbase.AbstractFourMomentum, mass::Real)
    return (mass * one(DiracMatrix) - slashed(mom)) / (sqrt(abs(getT(mom)) + mass))
end

function base_state(
    particle::Fermion,
    ::Incoming,
    mom::QEDbase.AbstractFourMomentum,
    spin::AbstractDefiniteSpin,
)
    booster = _booster_fermion(mom, mass(particle))
    return BiSpinor(booster[:, _spin_index(spin)])
end

function base_state(
    particle::Fermion, ::Incoming, mom::QEDbase.AbstractFourMomentum, spin::AllSpin
)
    booster = _booster_fermion(mom, mass(particle))
    return SVector(BiSpinor(booster[:, 1]), BiSpinor(booster[:, 2]))
end

function base_state(
    particle::AntiFermion,
    ::Incoming,
    mom::QEDbase.AbstractFourMomentum,
    spin::AbstractDefiniteSpin,
)
    booster = _booster_antifermion(mom, mass(particle))
    return AdjointBiSpinor(BiSpinor(booster[:, _spin_index(spin) + 2])) * GAMMA[1]
end

function base_state(
    particle::AntiFermion, ::Incoming, mom::QEDbase.AbstractFourMomentum, spin::AllSpin
)
    booster = _booster_antifermion(mom, mass(particle))
    return SVector(
        AdjointBiSpinor(BiSpinor(booster[:, 3])) * GAMMA[1],
        AdjointBiSpinor(BiSpinor(booster[:, 4])) * GAMMA[1],
    )
end

function base_state(
    particle::Fermion,
    ::Outgoing,
    mom::QEDbase.AbstractFourMomentum,
    spin::AbstractDefiniteSpin,
)
    booster = _booster_fermion(mom, mass(particle))
    return AdjointBiSpinor(BiSpinor(booster[:, _spin_index(spin)])) * GAMMA[1]
end

function base_state(
    particle::Fermion, ::Outgoing, mom::QEDbase.AbstractFourMomentum, spin::AllSpin
)
    booster = _booster_fermion(mom, mass(particle))
    return SVector(
        AdjointBiSpinor(BiSpinor(booster[:, 1])) * GAMMA[1],
        AdjointBiSpinor(BiSpinor(booster[:, 2])) * GAMMA[1],
    )
end

function base_state(
    particle::AntiFermion,
    ::Outgoing,
    mom::QEDbase.AbstractFourMomentum,
    spin::AbstractDefiniteSpin,
)
    booster = _booster_antifermion(mom, mass(particle))
    return BiSpinor(booster[:, _spin_index(spin) + 2])
end

function base_state(
    particle::AntiFermion, ::Outgoing, mom::QEDbase.AbstractFourMomentum, spin::AllSpin
)
    booster = _booster_antifermion(mom, mass(particle))
    return SVector(BiSpinor(booster[:, 3]), BiSpinor(booster[:, 4]))
end

function _photon_state(pol::AllPolarization, mom::QEDbase.AbstractFourMomentum)
    cth = getCosTheta(mom)
    sth = sqrt(1 - cth^2)
    cos_phi = getCosPhi(mom)
    sin_phi = getSinPhi(mom)
    return SVector(
        SLorentzVector{Float64}(0.0, cth * cos_phi, cth * sin_phi, -sth),
        SLorentzVector{Float64}(0.0, -sin_phi, cos_phi, 0.0),
    )
end

function _photon_state(pol::PolarizationX, mom::QEDbase.AbstractFourMomentum)
    cth = getCosTheta(mom)
    sth = sqrt(1 - cth^2)
    cos_phi = getCosPhi(mom)
    sin_phi = getSinPhi(mom)
    return SLorentzVector{Float64}(0.0, cth * cos_phi, cth * sin_phi, -sth)
end

function _photon_state(pol::PolarizationY, mom::QEDbase.AbstractFourMomentum)
    cos_phi = getCosPhi(mom)
    sin_phi = getSinPhi(mom)
    return SLorentzVector{Float64}(0.0, -sin_phi, cos_phi, 0.0)
end

@inline function base_state(
    particle::Photon,
    ::ParticleDirection,
    mom::QEDbase.AbstractFourMomentum,
    pol::AllPolarization,
)
    return _photon_state(pol, mom)
end

@inline function base_state(
    particle::Photon,
    ::ParticleDirection,
    mom::QEDbase.AbstractFourMomentum,
    pol::AbstractPolarization,
)
    return _photon_state(pol, mom)
end
