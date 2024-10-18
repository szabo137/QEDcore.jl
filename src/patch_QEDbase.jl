
#############
# Interface: Phase space layout
#############

export AbstractPhaseSpaceLayout
export AbstractInPhaseSpaceLayout
export AbstractOutPhaseSpaceLayout
export build_momenta

abstract type AbstractPhaseSpaceLayout end
abstract type AbstractInPhaseSpaceLayout <: AbstractPhaseSpaceLayout end
abstract type AbstractOutPhaseSpaceLayout{InPSL<:AbstractInPhaseSpaceLayout} <:
              AbstractPhaseSpaceLayout end

"""

    phase_space_dimension(proc,model,::AbstractPhaseSpaceLayout)

TBW
"""
function phase_space_dimension end

"""

    in_phase_space_layout(::AbstactOutPhaseSpaceLayout)

TBW
"""
function in_phase_space_layout end

"""

    _build_momenta(proc,model,in_psl, in_coords)
    _build_momenta(proc,model,Ptot,out_psl,out_coords)

TBW
"""
function _build_momenta end

#############
# Implementations: Phase space layout
#############

### in ps layout

function _build_momenta(
    ::Val{Nc},
    proc::AbstractProcessDefinition,
    model::AbstractModelDefinition,
    in_psl::AbstractInPhaseSpaceLayout,
    in_coords::NTuple{N},
) where {Nc,N}
    throw(
        InvalidInputError(
            "number of coordinates <$N> must be the same as the phase-space dimension <$Nc>"
        ),
    )
end

function _build_momenta(
    ::Val{Nc},
    proc::AbstractProcessDefinition,
    model::AbstractModelDefinition,
    in_psl::AbstractInPhaseSpaceLayout,
    in_coords::NTuple{Nc,T},
) where {Nc,T}
    return _build_momenta(proc, model, in_psl, in_coords)
end

function build_momenta(
    proc::AbstractProcessDefinition,
    model::AbstractModelDefinition,
    in_psl::AbstractInPhaseSpaceLayout,
    in_coords::Tuple,
)
    return _build_momenta(
        Val(phase_space_dimension(proc, model, in_psl)), proc, model, in_psl, in_coords
    )
end

"""
Scalar version of _build_momenta for in_psl
"""
function build_momenta(
    proc::AbstractProcessDefinition,
    model::AbstractModelDefinition,
    in_psl::AbstractInPhaseSpaceLayout,
    in_coords::Real,
)
    return build_momenta(proc, model, in_psl, (in_coords,))
end

### out ps layout

function _build_momenta(
    ::Val{Nc},
    proc::AbstractProcessDefinition,
    model::AbstractModelDefinition,
    Ptot::AbstractFourMomentum,
    out_psl::AbstractOutPhaseSpaceLayout,
    out_coords::NTuple{N},
) where {Nc,N}
    throw(
        InvalidInputError(
            "number of coordinates <$N> must be the same as the out-phase-space dimension <$Nc>",
        ),
    )
end

function _build_momenta(
    ::Val{Nc},
    proc::AbstractProcessDefinition,
    model::AbstractModelDefinition,
    Ptot::AbstractFourMomentum,
    out_psl::AbstractOutPhaseSpaceLayout,
    out_coords::NTuple{Nc,T},
) where {Nc,T<:Real}
    return _build_momenta(proc, model, Ptot, out_psl, out_coords)
end

function build_momenta(
    proc::AbstractProcessDefinition,
    model::AbstractModelDefinition,
    Ptot::AbstractFourMomentum,
    out_psl::AbstractOutPhaseSpaceLayout,
    out_coords::NTuple{Nc,T},
) where {Nc,T}
    return _build_momenta(
        Val(phase_space_dimension(proc, model, out_psl)),
        proc,
        model,
        Ptot,
        out_psl,
        out_coords,
    )
end
