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
SCREENINGS = (NoScreening(), Screening())


@testset "ne = $ne_ccm" for ne_ccm in NES_ccm
    ne_internal = MatterModels._internalize_density(ne_ccm)

    KF = MatterModels._fermi_wave_vector(ne_internal)
    EF = MatterModels._fermi_energy_from_kF(KF)
    N0 = KF / (2 * pi^2)

    OMS = EF .* (0.0, rand(RNG), 1 + rand(RNG), 2 + rand(RNG), 3 + rand(RNG))
    QS = KF .* (1.0e-2 * rand(RNG), rand(RNG), 2 + rand(RNG), 3 + rand(RNG))

    ### Zero temperature

    @testset "scr = $scr" for scr in SCREENINGS
        @testset "T = 0.0 eV" begin
            test_proper_system = IdealElectronSystem{ZeroTemperature}(ne_ccm)
            test_system = InteractingElectronSystem(test_proper_system, scr)

            @testset "properties" begin
                _rf_property_check(test_system, ne_internal, zero(ne_internal))
                @test screening(test_system) == scr
                @test iszero(temperature(test_system))
                @test isinf(beta(test_system))
                @test isinf(betabar(test_system))
            end

            @testset "q: $q" for q in QS

                @testset "rf stability" begin
                    _rf_stability_check(test_system, q)
                end

                # TODO: find out, why f-sum rule does not work for screened rf
                if scr == NoScreening()
                    @testset "rf f-sum rule" begin
                        _f_sum_rule_check(test_system, q)
                    end
                end

                @testset "om: $om" for om in OMS

                    @testset "rf symmetry" begin
                        _rf_symmetry_check(test_system, om, q, RTOL)
                    end

                    @testset "rf sanity check" begin
                        groundtruth_imag_rf = imag(dynamic_response(test_system, (om, q)))
                        @test isapprox(
                            groundtruth_imag_rf,
                            imag_dynamic_response(test_system, (om, q)),
                            rtol = RTOL,
                        )

                        groundtruth_real_rf = real(dynamic_response(test_system, (om, q)))

                        @test isapprox(
                            groundtruth_real_rf,
                            real_dynamic_response(test_system, (om, q)),
                            rtol = RTOL,
                        )

                        groundtruth_rf =
                            real_dynamic_response(test_system, (om, q)) +
                            1im * imag_dynamic_response(test_system, (om, q))
                        @test isapprox(
                            groundtruth_rf,
                            dynamic_response(test_system, (om, q)),
                            rtol = RTOL,
                        )

                        proper_rf = dynamic_response(test_proper_system, (om, q))
                        lep = local_effective_potential(test_system, (om, q))
                        groundtruth_rf_from_proper = proper_rf / (1 - lep * proper_rf)
                        @test isapprox(
                            groundtruth_rf_from_proper,
                            dynamic_response(test_system, (om, q)),
                            rtol = RTOL,
                        )

                        diel_func = dielectric_function(test_system, (om, q))
                        groundtruth_rf_from_dielec_func = proper_rf / diel_func
                        @test isapprox(
                            groundtruth_rf_from_dielec_func,
                            dynamic_response(test_system, (om, q)),
                            rtol = RTOL,
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

                test_proper_system = IdealElectronSystem(ne_ccm, T_eV, approx)
                test_system = InteractingElectronSystem(test_proper_system, scr)

                @testset "properties" begin
                    _rf_property_check(test_system, ne_internal, T_internal)
                    #@test response_approximation(test_system) == approx
                    #@test screening(test_system) == scr
                end

                @testset "q: $q" for q in QS

                    @testset "rf stability" begin
                        _rf_stability_check(test_system, q)
                    end

                    #=
                    if approx == NoApprox()
                        @testset "rf f-sum rule" begin
                            _f_sum_rule_check(test_system, q)
                        end
                    end
                    =#

                    @testset "om: $om" for om in OMS

                        @testset "rf symmetry" begin
                            _rf_symmetry_check(test_system, om, q, RTOL)
                        end

                        @testset "rf sanity check" begin
                            groundtruth_imag_rf = imag(dynamic_response(test_system, (om, q)))
                            @test isapprox(
                                groundtruth_imag_rf,
                                imag_dynamic_response(test_system, (om, q)),
                                rtol = RTOL,
                            )

                            groundtruth_real_rf = real(dynamic_response(test_system, (om, q)))

                            @test isapprox(
                                groundtruth_real_rf,
                                real_dynamic_response(test_system, (om, q)),
                                rtol = RTOL,
                            )

                            groundtruth_rf =
                                real_dynamic_response(test_system, (om, q)) +
                                1im * imag_dynamic_response(test_system, (om, q))
                            @test isapprox(
                                groundtruth_rf,
                                dynamic_response(test_system, (om, q)),
                                rtol = RTOL,
                            )

                            proper_rf = dynamic_response(test_proper_system, (om, q))
                            lep = local_effective_potential(test_system, (om, q))
                            groundtruth_rf_from_proper = proper_rf / (1 - lep * proper_rf)
                            @test isapprox(
                                groundtruth_rf_from_proper,
                                dynamic_response(test_system, (om, q)),
                                rtol = RTOL,
                            )

                            diel_func = dielectric_function(test_system, (om, q))
                            groundtruth_rf_from_dielec_func = proper_rf / diel_func
                            @test isapprox(
                                groundtruth_rf_from_dielec_func,
                                dynamic_response(test_system, (om, q)),
                                rtol = RTOL,
                            )
                        end

                        @testset "dsf sanity check" begin

                            # TODO: find out, why this is broken sometimes
                            #=
                            if !iszero(om)
                                _dsf_sanity_check(test_system, om, q)
                            end
                            =#
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
    end # screening
end
