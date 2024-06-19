using QEDbase: QEDbase
using QEDcore
using StaticArrays
using Random

include("../utils.jl")

FERMION_STATES_GROUNDTRUTH_FACTORY = Dict(
    (QEDbase.Incoming, Electron) => IncomingFermionSpinor,
    (QEDbase.Outgoing, Electron) => OutgoingFermionSpinor,
    (QEDbase.Incoming, Positron) => IncomingAntiFermionSpinor,
    (QEDbase.Outgoing, Positron) => OutgoingAntiFermionSpinor,
)

RNG = MersenneTwister(708583836976)
ATOL = eps()
RTOL = 0.0
PHOTON_ENERGIES = (0.0, rand(RNG), rand(RNG) * 10)
COS_THETAS = (-1.0, -rand(RNG), 0.0, rand(RNG), 1.0)

# check every quadrant
PHIS = (
    0.0,
    rand(RNG) * pi / 2,
    pi / 2,
    (1.0 + rand(RNG)) * pi / 2,
    pi,
    (2 + rand(RNG)) * pi / 2,
    3 * pi / 2,
    (3 + rand(RNG)) * pi / 2,
    2 * pi,
)

X, Y, Z = rand(RNG, 3)

# test function to test scalar broadcasting
test_broadcast(x::QEDbase.AbstractParticle) = x
test_broadcast(x::QEDbase.ParticleDirection) = x
test_broadcast(x::QEDbase.AbstractSpinOrPolarization) = x

@testset "fermion" begin
    struct TestFermion <: Fermion end
    @test QEDbase.is_fermion(TestFermion())
    @test QEDbase.is_particle(TestFermion())
    @test !QEDbase.is_anti_particle(TestFermion())
    @test test_broadcast.(TestFermion()) == TestFermion()

    @testset "$p $d" for (p, d) in Iterators.product(
        (Electron, Positron), (QEDbase.Incoming, QEDbase.Outgoing)
    )
        mom = SFourMomentum(sqrt(QEDbase.mass(p()) + X^2 + Y^2 + Z^2), X, Y, Z)
        particle_mass = QEDbase.mass(p())
        groundtruth_states = FERMION_STATES_GROUNDTRUTH_FACTORY[(d, p)](mom, particle_mass)
        groundtruth_tuple = SVector(groundtruth_states(1), groundtruth_states(2))
        @test base_state(p(), d(), mom, QEDbase.AllSpin()) == groundtruth_tuple
        @test base_state(p(), d(), mom, QEDbase.SpinUp()) == groundtruth_tuple[1]
        @test base_state(p(), d(), mom, QEDbase.SpinDown()) == groundtruth_tuple[2]

        @test QEDbase._as_svec(base_state(p(), d(), mom, QEDbase.AllSpin())) isa SVector
        @test QEDbase._as_svec(base_state(p(), d(), mom, QEDbase.SpinUp())) isa SVector
        @test QEDbase._as_svec(base_state(p(), d(), mom, QEDbase.SpinDown())) isa SVector

        @test QEDbase._as_svec(base_state(p(), d(), mom, QEDbase.AllSpin())) ==
            groundtruth_tuple
        @test QEDbase._as_svec(base_state(p(), d(), mom, QEDbase.SpinUp()))[1] ==
            groundtruth_tuple[1]
        @test QEDbase._as_svec(base_state(p(), d(), mom, QEDbase.SpinDown()))[1] ==
            groundtruth_tuple[2]
    end
end

@testset "photon" begin
    @test !QEDbase.is_fermion(Photon())
    @test QEDbase.is_boson(Photon())
    @test QEDbase.is_particle(Photon())
    @test QEDbase.is_anti_particle(Photon())
    @test QEDbase.charge(Photon()) == 0.0
    @test QEDbase.mass(Photon()) == 0.0
    @test test_broadcast.(Photon()) == Photon()

    @testset "$D" for D in [QEDbase.Incoming, QEDbase.Outgoing]
        @testset "$om $cth $phi" for (om, cth, phi) in
                                     Iterators.product(PHOTON_ENERGIES, COS_THETAS, PHIS)
            #@testset "$x $y $z" for (x,y,z) in Iterators.product(X_arr,Y_arr,Z_arr)

            mom = SFourMomentum(_cartesian_coordinates(om, om, cth, phi))
            both_photon_states = base_state(Photon(), D(), mom, QEDbase.AllPolarization())

            # property test the photon states
            @test isapprox((both_photon_states[1] * mom), 0.0, atol=ATOL, rtol=RTOL)
            @test isapprox((both_photon_states[2] * mom), 0.0, atol=ATOL, rtol=RTOL)
            @test isapprox(
                (both_photon_states[1] * both_photon_states[1]), -1.0, atol=ATOL, rtol=RTOL
            )
            @test isapprox(
                (both_photon_states[2] * both_photon_states[2]), -1.0, atol=ATOL, rtol=RTOL
            )
            @test isapprox(
                (both_photon_states[1] * both_photon_states[2]), 0.0, atol=ATOL, rtol=RTOL
            )

            # test the single polarization states
            @test base_state(Photon(), D(), mom, QEDbase.PolarizationX()) ==
                both_photon_states[1]
            @test base_state(Photon(), D(), mom, QEDbase.PolarizationY()) ==
                both_photon_states[2]
            @test base_state(Photon(), D(), mom, QEDbase.PolX()) == both_photon_states[1]
            @test base_state(Photon(), D(), mom, QEDbase.PolY()) == both_photon_states[2]

            @test QEDbase._as_svec(base_state(Photon(), D(), mom, QEDbase.PolX())) isa
                SVector
            @test QEDbase._as_svec(base_state(Photon(), D(), mom, QEDbase.PolY())) isa
                SVector
            @test QEDbase._as_svec(base_state(Photon(), D(), mom, QEDbase.AllPol())) isa
                SVector

            @test QEDbase._as_svec(base_state(Photon(), D(), mom, QEDbase.PolX()))[1] ==
                both_photon_states[1]
            @test QEDbase._as_svec(base_state(Photon(), D(), mom, QEDbase.PolY()))[1] ==
                both_photon_states[2]
            @test QEDbase._as_svec(base_state(Photon(), D(), mom, QEDbase.AllPol()))[1] ==
                both_photon_states[1]
            @test QEDbase._as_svec(base_state(Photon(), D(), mom, QEDbase.AllPol()))[2] ==
                both_photon_states[2]
        end
    end
end
