#######
#
# Abstract types
#
#######

import QEDbase: getT, getX, getY, getZ, setT!, setX!, setY!, setZ!
import StaticArrays: similar_type

#######
#
# Concrete LorentzVector types
#
#######
"""
$(TYPEDEF)

Concrete implementation of a generic static Lorentz vector. Each manipulation of an concrete implementation which is not self-contained (i.e. produces the same Lorentz vector type) will result in this type.

# Fields
$(TYPEDFIELDS)
"""
struct SLorentzVector{T} <: QEDbase.AbstractLorentzVector{T}
    "`t` component"
    t::T

    "`x` component"
    x::T

    "`y` component"
    y::T

    "`z` component"
    z::T
end
SLorentzVector(t, x, y, z) = SLorentzVector(promote(t, x, y, z)...)

function similar_type(::Type{A}, ::Type{T}, ::Size{S}) where {A<:SLorentzVector,T,S}
    return SLorentzVector{T}
end

@inline getT(lv::SLorentzVector) = lv.t
@inline getX(lv::SLorentzVector) = lv.x
@inline getY(lv::SLorentzVector) = lv.y
@inline getZ(lv::SLorentzVector) = lv.z

# TODO: this breaks incremental compilation because it's trying to eval permanent changes in a different module
#register_LorentzVectorLike(SLorentzVector)
@traitimpl QEDbase.IsLorentzVectorLike{SLorentzVector}

"""
$(TYPEDEF)

Concrete implementation of a generic mutable Lorentz vector. Each manipulation of an concrete implementation which is not self-contained (i.e. produces the same Lorentz vector type) will result in this type.

# Fields
$(TYPEDFIELDS)
"""
mutable struct MLorentzVector{T} <: QEDbase.AbstractLorentzVector{T}
    "`t` component"
    t::T

    "`x` component"
    x::T

    "`y` component"
    y::T

    "`z` component"
    z::T
end
MLorentzVector(t, x, y, z) = MLorentzVector(promote(t, x, y, z)...)

function similar_type(::Type{A}, ::Type{T}, ::Size{S}) where {A<:MLorentzVector,T,S}
    return MLorentzVector{T}
end

@inline getT(lv::MLorentzVector) = lv.t
@inline getX(lv::MLorentzVector) = lv.x
@inline getY(lv::MLorentzVector) = lv.y
@inline getZ(lv::MLorentzVector) = lv.z

function setT!(lv::MLorentzVector, value::T) where {T}
    return lv.t = value
end

function setX!(lv::MLorentzVector, value::T) where {T}
    return lv.x = value
end

function setY!(lv::MLorentzVector, value::T) where {T}
    return lv.y = value
end

function setZ!(lv::MLorentzVector, value::T) where {T}
    return lv.z = value
end

# TODO: this breaks incremental compilation because it's trying to eval permanent changes in a different module
#register_LorentzVectorLike(MLorentzVector)
@traitimpl QEDbase.IsLorentzVectorLike{MLorentzVector}
@traitimpl QEDbase.IsMutableLorentzVectorLike{MLorentzVector}
