module QEDcore

import Base: *
import StaticArrays: similar_type

# lorentz vectors
export SLorentzVector, MLorentzVector

# four momenta
export SFourMomentum, MFourMomentum

# spinors
export BiSpinor, AdjointBiSpinor, DiracMatrix

# gamma matrices
export gamma, GAMMA, DiracGammaRepresentation, slashed

# particle types
export AbstractParticleSpinor
export FermionLike, Fermion, AntiFermion, MajoranaFermion
export BosonLike, Boson, AntiBoson, MajoranaBoson
export Electron, Positron, Photon

# particle spinors
export BASE_PARTICLE_SPINOR, BASE_ANTIPARTICLE_SPINOR
export IncomingFermionSpinor,
    OutgoingFermionSpinor, IncomingAntiFermionSpinor, OutgoingAntiFermionSpinor
export SpinorU, SpinorUbar, SpinorV, SpinorVbar
export @valid_spinor_input

using QEDbase: QEDbase
using DocStringExtensions
using StaticArrays
using SimpleTraits

include("algebraic_objects/dirac_tensors/types.jl")
include("algebraic_objects/dirac_tensors/multiplication.jl")

include("algebraic_objects/four_momentum.jl")
include("algebraic_objects/lorentz_vector.jl")
include("algebraic_objects/gamma_matrices.jl")

include("particles/particle_types.jl")
include("particles/spinors.jl")

end
