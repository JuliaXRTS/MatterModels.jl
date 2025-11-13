function imag_dynamic_response(
        elsys::AbstractIdealElectronSystem,
        om_q::NTuple{2, T},
    ) where {T <: Real}
    ombar, qbar = _transform_om_q(elsys, om_q)
    approx = response_approximation(elsys)
    KF = fermi_wave_vector(elsys)
    N0 = KF / (2 * pi^2)


    if ombar <= zero(ombar)
        if qbar < zero(qbar)
            return N0 * _imag_ideal_dynamic_response(elsys, approx, -ombar, -qbar)
        else
            return -N0 * _imag_ideal_dynamic_response(elsys, approx, -ombar, qbar)
        end
    end

    return N0 * _imag_ideal_dynamic_response(elsys, approx, ombar, qbar)
end

function real_dynamic_response(
        elsys::AbstractIdealElectronSystem,
        om_q::NTuple{2, T},
    ) where {T <: Real}
    ombar, qbar = _transform_om_q(elsys, om_q)

    # TODO: insert symmetries here as well

    KF = fermi_wave_vector(elsys)
    N0 = KF / (2 * pi^2)

    if ombar <= zero(ombar)
        return N0 * _real_ideal_dynamic_response(elsys, response_approximation(elsys), -ombar, qbar)
    end

    return N0 * _real_ideal_dynamic_response(elsys, response_approximation(elsys), ombar, qbar)
end
