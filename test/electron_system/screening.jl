using Test
using Random
using QuadGK
using Unitful

using MatterModels

include("checks.jl")

RNG = Xoshiro(137)
ATOL = sqrt(eps())
RTOL = 1.0e-4

_transform12(x01) = x01 + one(x01)

NES_ccm = 1u"cm^(-3)" .* (
    _transform12(rand(RNG)) * 1.0e20,
    _transform12(rand(RNG)) * 1.0e21,
    _transform12(rand(RNG)) * 1.0e22,
    _transform12(rand(RNG)) * 1.0e23,
    _transform12(rand(RNG)) * 1.0e24,
)
TEMPS_eV = 1u"eV" .* [
    rand(RNG),
    rand(RNG) * 10,
    rand(RNG) * 100,
    rand(RNG) * 1000,
]

APPROXS = [NoApprox(), NonDegenerated(), Degenerated()]


@testset "ne = $ne_ccm" for ne_ccm in NES_ccm
    ne_internal = MatterModels._internalize_density(ne_ccm)

    KF = MatterModels._fermi_wave_vector(ne_internal)
    EF = MatterModels._fermi_energy_from_kF(KF)
    N0 = KF / (2 * pi^2)

    OMS = EF .* (0.0, rand(RNG), 1 + rand(RNG), 2 + rand(RNG), 3 + rand(RNG))
    QS = KF .* (1.0e-2 * rand(RNG), rand(RNG), 2 + rand(RNG), 3 + rand(RNG))

    TEST_SYSTEM_zeroT = IdealElectronSystem{ZeroTemperature}(ne_ccm)
    TEST_SYSTEMS_finT = [IdealElectronSystem(ne_ccm, T_eV, approx) for T_eV in TEMPS_eV for approx in APPROXS]
    TEST_SYSTEMS = (TEST_SYSTEM_zeroT, TEST_SYSTEMS_finT...)

    @testset "T = $(temperature(test_system))" for test_system in TEST_SYSTEMS
        @testset "no screening" begin
            scr = NoScreening()
            @testset "om: $om, q: $q" for (om, q) in Iterators.product(OMS, QS)

                @test isapprox(dielectric_function(test_system, scr, (om, q)), one(om))

                @test isapprox(local_effective_potential(scr, (om, q)), zero(om))

                @test isapprox(pseudo_potential(scr, (om, q)), zero(om))

                @test isapprox(local_field_correction(scr, (om, q)), zero(om))
            end
        end

        @testset "screening" begin
            # checking only default value: coulomb+no_lfc
            scr = Screening()

            @testset "default constructor" begin
                @test scr == Screening(CoulombPseudoPotential(), NoLocalFieldCorrection())
            end

            @testset "om: $om, q: $q" for (om, q) in Iterators.product(OMS, QS)

                @test isapprox(
                    dielectric_function(test_system, scr, (om, q)),
                    one(om) -
                        pseudo_potential(scr, (om, q)) * dynamic_response(test_system, (om, q)),
                )

                @test isapprox(
                    local_effective_potential(scr, (om, q)),
                    pseudo_potential(scr, (om, q)),
                )

                @test isapprox(pseudo_potential(scr, (om, q)), ELEMENTARY_CHARGE_SQUARED / q^2)

                @test isapprox(local_field_correction(scr, (om, q)), zero(om))
            end
        end

    end

end
