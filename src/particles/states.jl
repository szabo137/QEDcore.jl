import QEDbase: base_state

function _booster_fermion(mom::QEDbase.AbstractFourMomentum, mass::Real)
    return (slashed(mom) + mass * one(DiracMatrix)) / (sqrt(abs(getT(mom)) + mass))
end

function _booster_antifermion(mom::QEDbase.AbstractFourMomentum, mass::Real)
    return (mass * one(DiracMatrix) - slashed(mom)) / (sqrt(abs(getT(mom)) + mass))
end

function QEDbase.base_state(
    particle::QEDbase.Fermion,
    ::QEDbase.Incoming,
    mom::QEDbase.AbstractFourMomentum,
    spin::QEDbase.AbstractDefiniteSpin,
)
    booster = _booster_fermion(mom, QEDbase.mass(particle))
    return BiSpinor(booster[:, QEDbase._spin_index(spin)])
end

function QEDbase.base_state(
    particle::QEDbase.Fermion,
    ::QEDbase.Incoming,
    mom::QEDbase.AbstractFourMomentum,
    spin::QEDbase.AllSpin,
)
    booster = _booster_fermion(mom, QEDbase.mass(particle))
    return SVector(BiSpinor(booster[:, 1]), BiSpinor(booster[:, 2]))
end

function QEDbase.base_state(
    particle::QEDbase.AntiFermion,
    ::QEDbase.Incoming,
    mom::QEDbase.AbstractFourMomentum,
    spin::QEDbase.AbstractDefiniteSpin,
)
    booster = _booster_antifermion(mom, QEDbase.mass(particle))
    return AdjointBiSpinor(BiSpinor(booster[:, QEDbase._spin_index(spin) + 2])) * GAMMA[1]
end

function QEDbase.base_state(
    particle::QEDbase.AntiFermion,
    ::QEDbase.Incoming,
    mom::QEDbase.AbstractFourMomentum,
    spin::QEDbase.AllSpin,
)
    booster = _booster_antifermion(mom, QEDbase.mass(particle))
    return SVector(
        AdjointBiSpinor(BiSpinor(booster[:, 3])) * GAMMA[1],
        AdjointBiSpinor(BiSpinor(booster[:, 4])) * GAMMA[1],
    )
end

function QEDbase.base_state(
    particle::QEDbase.Fermion,
    ::QEDbase.Outgoing,
    mom::QEDbase.AbstractFourMomentum,
    spin::QEDbase.AbstractDefiniteSpin,
)
    booster = _booster_fermion(mom, QEDbase.mass(particle))
    return AdjointBiSpinor(BiSpinor(booster[:, QEDbase._spin_index(spin)])) * GAMMA[1]
end

function QEDbase.base_state(
    particle::QEDbase.Fermion,
    ::QEDbase.Outgoing,
    mom::QEDbase.AbstractFourMomentum,
    spin::QEDbase.AllSpin,
)
    booster = _booster_fermion(mom, QEDbase.mass(particle))
    return SVector(
        AdjointBiSpinor(BiSpinor(booster[:, 1])) * GAMMA[1],
        AdjointBiSpinor(BiSpinor(booster[:, 2])) * GAMMA[1],
    )
end

function QEDbase.base_state(
    particle::QEDbase.AntiFermion,
    ::QEDbase.Outgoing,
    mom::QEDbase.AbstractFourMomentum,
    spin::QEDbase.AbstractDefiniteSpin,
)
    booster = _booster_antifermion(mom, QEDbase.mass(particle))
    return BiSpinor(booster[:, QEDbase._spin_index(spin) + 2])
end

function QEDbase.base_state(
    particle::QEDbase.AntiFermion,
    ::QEDbase.Outgoing,
    mom::QEDbase.AbstractFourMomentum,
    spin::QEDbase.AllSpin,
)
    booster = _booster_antifermion(mom, QEDbase.mass(particle))
    return SVector(BiSpinor(booster[:, 3]), BiSpinor(booster[:, 4]))
end

function _photon_state(pol::QEDbase.AllPolarization, mom::QEDbase.AbstractFourMomentum)
    cth = QEDbase.getCosTheta(mom)
    sth = sqrt(1 - cth^2)
    cos_phi = QEDbase.getCosPhi(mom)
    sin_phi = QEDbase.getSinPhi(mom)
    return SVector(
        SLorentzVector{Float64}(0.0, cth * cos_phi, cth * sin_phi, -sth),
        SLorentzVector{Float64}(0.0, -sin_phi, cos_phi, 0.0),
    )
end

function _photon_state(pol::QEDbase.PolarizationX, mom::QEDbase.AbstractFourMomentum)
    cth = QEDbase.getCosTheta(mom)
    sth = sqrt(1 - cth^2)
    cos_phi = QEDbase.getCosPhi(mom)
    sin_phi = QEDbase.getSinPhi(mom)
    return SLorentzVector{Float64}(0.0, cth * cos_phi, cth * sin_phi, -sth)
end

function _photon_state(pol::QEDbase.PolarizationY, mom::QEDbase.AbstractFourMomentum)
    cos_phi = QEDbase.getCosPhi(mom)
    sin_phi = QEDbase.getSinPhi(mom)
    return SLorentzVector{Float64}(0.0, -sin_phi, cos_phi, 0.0)
end

@inline function QEDbase.base_state(
    particle::QEDbase.Photon,
    ::QEDbase.ParticleDirection,
    mom::QEDbase.AbstractFourMomentum,
    pol::QEDbase.AllPolarization,
)
    return _photon_state(pol, mom)
end

@inline function QEDbase.base_state(
    particle::QEDbase.Photon,
    ::QEDbase.ParticleDirection,
    mom::QEDbase.AbstractFourMomentum,
    pol::QEDbase.AbstractPolarization,
)
    return _photon_state(pol, mom)
end
