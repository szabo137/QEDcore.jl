#######
#
# Concrete Dirac Tensor types
#
#######
"""
$(TYPEDEF)

Concrete type to model a Dirac four-spinor with complex-valued components. These are the elements of an actual spinor space.
"""
struct BiSpinor <: QEDbase.AbstractDiracVector{ComplexF64}
    el1::ComplexF64
    el2::ComplexF64
    el3::ComplexF64
    el4::ComplexF64
end

"""
$(TYPEDEF)

Concrete type to model an adjoint Dirac four-spinor with complex-valued components. These are the elements of the dual spinor space.
"""
struct AdjointBiSpinor <: QEDbase.AbstractDiracVector{ComplexF64}
    el1::ComplexF64
    el2::ComplexF64
    el3::ComplexF64
    el4::ComplexF64
end

#interface
AdjointBiSpinor(spn::BiSpinor) = AdjointBiSpinor(conj(SVector(spn)))
BiSpinor(spn::AdjointBiSpinor) = BiSpinor(conj(SVector(spn)))

"""
$(TYPEDEF)

Concrete type to model Dirac matrices, i.e. matrix representations of linear mappings between two spinor spaces.
"""
struct DiracMatrix <: QEDbase.AbstractDiracMatrix{ComplexF64}
    el11::ComplexF64
    el12::ComplexF64
    el13::ComplexF64
    el14::ComplexF64
    el21::ComplexF64
    el22::ComplexF64
    el23::ComplexF64
    el24::ComplexF64
    el31::ComplexF64
    el32::ComplexF64
    el33::ComplexF64
    el34::ComplexF64
    el41::ComplexF64
    el42::ComplexF64
    el43::ComplexF64
    el44::ComplexF64
end
