#######
#
# Concrete SFourMomentum type
#
#######

import QEDbase: getT, getX, getY, getZ, setT!, setX!, setY!, setZ!
import StaticArrays: similar_type

"""
$(TYPEDEF)

Builds a static LorentzVectorLike with real components used to statically model the four-momentum of a particle or field.

# Fields
$(TYPEDFIELDS)
"""
struct SFourMomentum <: QEDbase.AbstractFourMomentum
    "energy component"
    E::Float64

    "`x` component"
    px::Float64

    "`y` component"
    py::Float64

    "`z` component"
    pz::Float64
end

"""
$(SIGNATURES)

The interface transforms each number-like input to float64:

$(TYPEDSIGNATURES)
"""
function SFourMomentum(t::T, x::T, y::T, z::T) where {T<:Union{Integer,Rational,Irrational}}
    return SFourMomentum(float(t), x, y, z)
end

function similar_type(::Type{A}, ::Type{T}, ::Size{S}) where {A<:SFourMomentum,T<:Real,S}
    return SFourMomentum
end
function similar_type(::Type{A}, ::Type{T}, ::Size{S}) where {A<:SFourMomentum,T,S}
    return SLorentzVector{T}
end

@inline getT(p::SFourMomentum) = p.E
@inline getX(p::SFourMomentum) = p.px
@inline getY(p::SFourMomentum) = p.py
@inline getZ(p::SFourMomentum) = p.pz

# TODO: this breaks incremental compilation because it's trying to eval permanent changes in a different module
#register_LorentzVectorLike(SFourMomentum)
@traitimpl QEDbase.IsLorentzVectorLike{SFourMomentum}

#######
#
# Concrete MFourMomentum type
#
#######
"""
$(TYPEDEF)

Builds a mutable LorentzVector with real components used to statically model the four-momentum of a particle or field.

# Fields
$(TYPEDFIELDS)
"""
mutable struct MFourMomentum <: QEDbase.AbstractFourMomentum
    "energy component"
    E::Float64

    "`x` component"
    px::Float64

    "`y` component"
    py::Float64

    "`z` component"
    pz::Float64
end

"""
$(SIGNATURES)

The interface transforms each number-like input to float64:

$(TYPEDSIGNATURES)
"""
function MFourMomentum(t::T, x::T, y::T, z::T) where {T<:Union{Integer,Rational,Irrational}}
    return MFourMomentum(float(t), x, y, z)
end

function similar_type(::Type{A}, ::Type{T}, ::Size{S}) where {A<:MFourMomentum,T<:Real,S}
    return MFourMomentum
end
function similar_type(::Type{A}, ::Type{T}, ::Size{S}) where {A<:MFourMomentum,T,S}
    return MLorentzVector{T}
end

@inline getT(p::MFourMomentum) = p.E
@inline getX(p::MFourMomentum) = p.px
@inline getY(p::MFourMomentum) = p.py
@inline getZ(p::MFourMomentum) = p.pz

function QEDbase.setT!(lv::MFourMomentum, value::Float64)
    return lv.E = value
end

function QEDbase.setX!(lv::MFourMomentum, value::Float64)
    return lv.px = value
end

function QEDbase.setY!(lv::MFourMomentum, value::Float64)
    return lv.py = value
end

function QEDbase.setZ!(lv::MFourMomentum, value::Float64)
    return lv.pz = value
end

# TODO: this breaks incremental compilation because it's trying to eval permanent changes in a different module
# register_LorentzVectorLike(MFourMomentum)
@traitimpl QEDbase.IsLorentzVectorLike{MFourMomentum}
@traitimpl QEDbase.IsMutableLorentzVectorLike{MFourMomentum}
