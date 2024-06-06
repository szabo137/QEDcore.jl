module QEDcore

import Base: *

# lorentz vectors
export SLorentzVector, MLorentzVector

# four momenta
export SFourMomentum, MFourMomentum

# spinors
export BiSpinor, AdjointBiSpinor, DiracMatrix

# gamma matrices
export gamma, GAMMA, DiracGammaRepresentation, slashed

# particle spinors
export BASE_PARTICLE_SPINOR, BASE_ANTIPARTICLE_SPINOR
export IncomingFermionSpinor,
    OutgoingFermionSpinor, IncomingAntiFermionSpinor, OutgoingAntiFermionSpinor
export SpinorU, SpinorUbar, SpinorV, SpinorVbar
export @valid_spinor_input

# phase space
export SphericalCoordinateSystem
export CenterOfMomentumFrame, ElectronRestFrame
export PhasespaceDefinition
export ParticleStateful, PhaseSpacePoint, InPhaseSpacePoint, OutPhaseSpacePoint
export spin, polarization, particle_direction, particle_species, momentum, momenta, getindex

using QEDbase
using DocStringExtensions
using StaticArrays
using SimpleTraits

include("dirac_tensors/types.jl")
include("dirac_tensors/multiplication.jl")

include("phase_spaces/types.jl")
include("phase_spaces/access.jl")
include("phase_spaces/create.jl")
include("phase_spaces/print.jl")
include("phase_spaces/utility.jl")

include("four_momentum.jl")

include("lorentz_vector.jl")

include("gamma_matrices.jl")

include("particles/states.jl")
include("particles/spinors.jl")

end
