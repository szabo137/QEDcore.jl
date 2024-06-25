import QEDbase:
    particle_direction,
    particle_species,
    momentum,
    process,
    model,
    phase_space_definition,
    particles

# accessor interface particle stateful
particle_direction(part::ParticleStateful) = part.dir
particle_species(part::ParticleStateful) = part.species
momentum(part::ParticleStateful) = part.mom

# accessor interface phase space point
"""
    Base.getindex(psp::PhaseSpacePoint, dir::Incoming, n::Int)

Overload for the array indexing operator `[]`. Returns the nth incoming particle in this phase space point.
"""
function Base.getindex(psp::PhaseSpacePoint, ::Incoming, n::Int)
    return psp.in_particles[n]
end

"""
    Base.getindex(psp::PhaseSpacePoint, dir::Outgoing, n::Int)

Overload for the array indexing operator `[]`. Returns the nth outgoing particle in this phase space point.
"""
function Base.getindex(psp::PhaseSpacePoint, ::Outgoing, n::Int)
    return psp.out_particles[n]
end

process(psp::PhaseSpacePoint) = psp.proc
model(psp::PhaseSpacePoint) = psp.model
phase_space_definition(psp::PhaseSpacePoint) = psp.ps_def

particles(psp::PhaseSpacePoint, ::Incoming) = psp.in_particles
particles(psp::PhaseSpacePoint, ::Outgoing) = psp.out_particles
