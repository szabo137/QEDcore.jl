### Trivial phase-space layouts

# maps all components onto four momenta
struct TrivialInPSL <: AbstractInPhaseSpaceLayout end

@inline QEDcore.phase_space_dimension(
    proc::AbstractProcessDefinition, ::AbstractModelDefinition, ::TrivialInPSL
) = 4 * number_incoming_particles(proc)

@inline function QEDcore._build_momenta(
    ::TestProcess, ::TestModel, ::TrivialInPSL, in_coords
)
    return _groundtruth_in_moms(in_coords)
end

# maps componets of N-1 particles onto four-momenta and uses energy-momentum conservation
struct TrivialOutPSL <: AbstractOutPhaseSpaceLayout{TrivialInPSL}
    in_psl::TrivialInPSL
end

@inline QEDcore.in_phase_space_layout(psl::TrivialOutPSL) = psl.in_psl

@inline QEDcore.phase_space_dimension(
    proc::AbstractProcessDefinition, ::AbstractModelDefinition, ::TrivialOutPSL
) = 4 * number_outgoing_particles(proc) - 4

@inline function QEDcore._build_momenta(
    proc::TestProcess,
    model::TestModel,
    Ptot::AbstractFourMomentum,
    out_psl::TrivialOutPSL,
    out_coords,
)
    return _groundtruth_out_moms(Ptot, out_coords)
end
