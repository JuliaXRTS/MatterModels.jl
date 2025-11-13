# checks for dynamic structure factors

function _dsf_detailed_balance_check(sys, om, q)

    input_pos_om = (om, q)
    input_neg_om = (-om, q)

    @test isapprox(
        dynamic_structure_factor(sys, input_neg_om),
        exp(-om * beta(sys)) *
            dynamic_structure_factor(sys, input_pos_om),
    )
    return nothing
end

function _dsf_sanity_check(sys, om, q)
    fac = pi * electron_density(sys)

    im_rf = imag_dynamic_response(sys, (om, q))
    T = temperature(sys)
    groundtruth_dsf =
        iszero(T) ? -im_rf / fac :
        -inv(fac * (one(om) - exp(-om / T))) * im_rf

    dsf = dynamic_structure_factor(sys, (om, q))
    @test isapprox(dsf, groundtruth_dsf)

    if dsf < 0.0
        @show dsf
        @show im_rf
    end
    @test dsf >= zero(dsf)

    return nothing
end

# checks for dynamic response functions

function _f_sum_rule_check(sys, q; bounded = true)
    T = temperature(sys)
    KF = fermi_wave_vector(sys)
    EF = fermi_energy(sys)
    N0 = KF / (2 * pi^2)

    qb = q / KF
    lower_om_bound = iszero(T) && bounded ? EF * max(zero(qb), qb^2 - 2 * qb) : zero(T)
    upper_om_bound = iszero(T) && bounded ? EF * (qb^2 + 2 * qb) : EF * (qb^2 + 2 * qb) * 300.0
    tmp, err = quadgk(
        x -> x * imag_dynamic_response(sys, (x, q)),
        lower_om_bound,
        upper_om_bound,
    )
    _first_moment = -2 * tmp / pi
    @test isapprox(
        _first_moment,
        electron_density(sys) * q^2,
        #4 * N0 * q^2 * EF / 3,
        rtol = 1.0e-2,
    )
    return nothing
end

function _rf_symmetry_check(sys, om, q, rtol)
    input_pos_q = (om, q)
    input_pos_om = (om, q)
    input_neg_om = (-om, q)
    input_neg_both = (-om, -q)

    #=
    # TODO: check if this is valid!

    @test isapprox(
        dynamic_response(sys, input_pos_q),
        dynamic_response(sys, input_neg_both),
        rtol = rtol,
    )
    =#
    @test isapprox(
        real_dynamic_response(sys, input_pos_om),
        real_dynamic_response(sys, input_neg_om),
        rtol = rtol,
    )
    @test isapprox(
        imag_dynamic_response(sys, input_pos_om),
        -imag_dynamic_response(sys, input_neg_om),
        rtol = rtol,
    )

    return nothing
end

function _rf_stability_check(sys, q)
    @test real_dynamic_response(sys, (0.0, q)) <= zero(q)
    @test imag_dynamic_response(sys, (0.0, q)) <= zero(q)
    return nothing
end

function _rf_property_check(sys, ne_internal, T_internal)
    KF_internal = cbrt(3 * pi^2 * ne_internal)
    EF_internal = KF_internal^2 / 2
    BETA_internal = inv(T_internal)
    BETABAR_internal = EF_internal * BETA_internal

    @test isapprox(temperature(sys), T_internal)
    @test isapprox(electron_density(sys), ne_internal)
    @test isapprox(beta(sys), BETA_internal)
    @test isapprox(betabar(sys), BETABAR_internal)
    @test isapprox(fermi_wave_vector(sys), KF_internal)
    @test isapprox(fermi_energy(sys), EF_internal)
    return nothing
end
