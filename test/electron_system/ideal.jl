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
TEMPS_eV = 1u"eV" .* (
    rand(RNG),
    rand(RNG) * 10,
    rand(RNG) * 100,
    rand(RNG) * 1000,
)

APPROXS = (NoApprox(), NonDegenerated(), Degenerated())


@testset "ne = $ne_ccm" for ne_ccm in NES_ccm
    ne_internal = MatterModels._internalize_density(ne_ccm)

    KF = MatterModels._fermi_wave_vector(ne_internal)
    EF = MatterModels._fermi_energy_from_kF(KF)
    N0 = KF / (2 * pi^2)

    OMS = EF .* (0.0, rand(RNG), 1 + rand(RNG), 2 + rand(RNG), 3 + rand(RNG))
    QS = KF .* (1.0e-2 * rand(RNG), rand(RNG), 2 + rand(RNG), 3 + rand(RNG))

    ### Zero temperature

    @testset "T = 0.0 eV" begin
        test_system = IdealElectronSystem{ZeroTemperature}(ne_ccm)
        test_system_small_finT = IdealElectronSystem(ne_ccm, 1.0e-2 * eps())

        @testset "properties" begin
            _rf_property_check(test_system, ne_internal, zero(ne_internal))
            @test iszero(temperature(test_system))
            @test isinf(beta(test_system))
            @test isinf(betabar(test_system))
        end

        @testset "q: $q" for q in QS

            @testset "rf stability" begin
                _rf_stability_check(test_system, q)
            end
            @testset "rf f-sum rule" begin
                _f_sum_rule_check(test_system, q)
            end

            @testset "om: $om" for om in OMS

                @testset "rf symmetry" begin
                    _rf_symmetry_check(test_system, om, q, RTOL)
                end

                @testset "rf sanity check" begin
                    groundtruth_imag_rf =
                        N0 * MatterModels._imag_lindhard_zero_temperature(NoApprox(), om / EF, q / KF)
                    @test isapprox(
                        groundtruth_imag_rf,
                        imag_dynamic_response(test_system, (om, q)),
                        rtol = RTOL,
                    )

                    groundtruth_real_rf =
                        N0 * MatterModels._real_lindhard_zero_temperature(NoApprox(), om / EF, q / KF)
                    @test isapprox(
                        groundtruth_real_rf,
                        real_dynamic_response(test_system, (om, q)),
                        rtol = RTOL,
                    )
                    groundtruth_rf = groundtruth_real_rf + 1im * groundtruth_imag_rf
                    @test isapprox(
                        groundtruth_rf,
                        dynamic_response(test_system, (om, q)),
                        rtol = RTOL,
                    )
                end

                @testset "rf matching" begin
                    @test isapprox(
                        imag_dynamic_response(test_system, (om, q)),
                        imag_dynamic_response(test_system_small_finT, (om, q)),
                        rtol = 1.0e-2,
                        atol = 1.0e-6,
                    )
                end

                @testset "dsf sanity check" begin
                    if !iszero(om)
                        _dsf_sanity_check(test_system, om, q)
                    end
                end
            end
        end
    end

    ### Finite temperature

    @testset "T = $T_eV" for T_eV in TEMPS_eV
        T_internal = MatterModels._internalize_temperature.(T_eV)
        @testset "approx = $approx" for approx in APPROXS

            test_system = IdealElectronSystem(ne_ccm, T_eV, approx)

            @testset "properties" begin
                _rf_property_check(test_system, ne_internal, T_internal)
                @test response_approximation(test_system) == approx
            end

            @testset "q: $q" for q in QS

                @testset "rf stability" begin
                    _rf_stability_check(test_system, q)
                end

                if approx == NoApprox()
                    @testset "rf f-sum rule" begin
                        _f_sum_rule_check(test_system, q)
                    end
                end

                @testset "om: $om" for om in OMS

                    @testset "rf symmetry" begin
                        _rf_symmetry_check(test_system, om, q, RTOL)
                    end

                    @testset "rf sanity check" begin
                        groundtruth_imag_rf = N0 * MatterModels._imag_lindhard_nonzero_temperature(
                            approx,
                            om / EF,
                            q / KF,
                            betabar(test_system),
                        )
                        @test isapprox(
                            groundtruth_imag_rf,
                            imag_dynamic_response(test_system, (om, q)),
                            rtol = RTOL,
                        )

                        groundtruth_real_rf = N0 * MatterModels._real_lindhard_nonzero_temperature(
                            approx,
                            om / EF,
                            q / KF,
                            betabar(test_system),
                        )
                        @test isapprox(
                            groundtruth_real_rf,
                            real_dynamic_response(test_system, (om, q)),
                        )

                        groundtruth_rf = groundtruth_real_rf + 1im * groundtruth_imag_rf
                        @test isapprox(groundtruth_rf, dynamic_response(test_system, (om, q)))
                    end

                    @testset "dsf sanity check" begin
                        if !iszero(om)
                            _dsf_sanity_check(test_system, om, q)
                        end
                    end

                    @testset "detailed balance" begin
                        if !iszero(om)
                            _dsf_detailed_balance_check(test_system, om, q)
                        end
                    end
                end # om
            end # q
        end # approx
    end # finite T

end
