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

function _imag_lindhard_zero_temperature(::NoApprox, ombar::T, qbar::T) where {T <: Real}
    return -pi / (2 * qbar) *
        (_imag_lindhardDSF_zeroT_minus(ombar, qbar) - _imag_lindhardDSF_zeroT_plus(ombar, qbar))
end

function _real_lindhard_zero_temperature(::NoApprox, ombar::T, qbar::T) where {T <: Real}
    num = _nu_minus(ombar, qbar)
    nup = _nu_plus(ombar, qbar)

    term1 = (1 - num^2) / 4 * _stable_log_term(num)
    term2 = (1 - nup^2) / 4 * _stable_log_term(nup)

    return -2 * (qbar / 2 - term1 + term2) / qbar
end

# finite temperature Lindhard, without any further approximation

## Real part

function _general_fermi(x, beta)
    mu = _chemical_potential_normalized(beta)
    denom = exp(beta * (x^2 - mu)) + one(beta)
    return inv(denom)
end

function _general_integrand(x, nu, beta)
    return x * _general_fermi(x, beta) * log(abs(x - nu) / abs(x + nu))
end

function _general_integal(nu, beta)
    res, _ = quadgk(x -> _general_integrand(x, nu, beta), 0.0, nu, Inf)

    return res
end

function _real_lindhard_nonzero_temperature(::NoApprox, ombar, qbar, bbar)
    num = _nu_minus(ombar, qbar)
    nup = _nu_plus(ombar, qbar)

    prefac = inv(qbar)

    res = -prefac * (_general_integal(num, bbar) - _general_integal(nup, bbar))

    return res
end

## Imaginary part

function _imag_lindhard_nonzero_temperature(::NoApprox, ombar::T, qbar::T, bbar::T) where {T <: Real}

    prefac = -pi / (2 * qbar)

    num = _nu_minus(ombar, qbar)
    nup = _nu_plus(ombar, qbar)
    mubar = _chemical_potential_normalized(bbar)

    term1 = log1pexp(bbar * (mubar - num^2))
    term2 = log1pexp(bbar * (mubar - nup^2))

    return prefac * (term1 - term2) / bbar

end
