# zero temperature Lindhard, without any further approximation

function _imag_lindhardDSF_zeroT_minus(omb::T, qb::T) where {T <: Real}
    num = _nu_minus(omb, qb)

    if abs(num) >= one(T)
        return zero(T)
    end
    return one(omb) - num^2

end

function _imag_lindhardDSF_zeroT_plus(omb::T, qb::T) where {T <: Real}
    nup = _nu_plus(omb, qb)
    if abs(nup) >= one(T)
        return zero(T)
    end
    return one(omb) - nup^2

end

function _imag_ideal_dynamic_response(::ZeroTemperatureApprox, ombar::T, qbar::T, betabar::T) where {T <: Real}
    return -pi / (2 * qbar) *
        (_imag_lindhardDSF_zeroT_minus(ombar, qbar) - _imag_lindhardDSF_zeroT_plus(ombar, qbar))
end

function _real_ideal_dynamic_response(::ZeroTemperatureApprox, ombar::T, qbar::T, betabar::T) where {T <: Real}
    num = _nu_minus(ombar, qbar)
    nup = _nu_plus(ombar, qbar)

    term1 = (1 - num^2) / 4 * _stable_log_term(num)
    term2 = (1 - nup^2) / 4 * _stable_log_term(nup)

    return -2 * (qbar / 2 - term1 + term2) / qbar
end
