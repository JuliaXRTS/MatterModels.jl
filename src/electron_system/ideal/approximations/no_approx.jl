# finite temperature Lindhard, without any further approximation

## Real part
# integraion method adopted from Eur. Phys. J. Plus (2016) 131: 114

function _fermi_prime(x, A, B)
    return term =
        return -2 * A * x / (
        exp(x^2 * A - B) * (exp(-x^2 * A + B) + 1)^2
    )
end

function _general_integrand(x, A, B)
    fp = _fermi_prime(x, A, B)
    return fp * (
        -x + (1 - x^2) * log(abs(x + 1) / abs(x - 1)) / 2
    )
end

function _general_integral(nu, beta)
    mu = ElectronicStructureModels._chemical_potential_normalized(beta)
    A = nu^2 * beta
    B = beta * mu
    res, _ = quadgk(x -> _general_integrand(x, A, B), 0.0, 1.0, Inf)

    return -nu^2 * res
end

function _real_ideal_dynamic_response(::NoApprox, ombar, qbar, bbar)
    num = _nu_minus(ombar, qbar)
    nup = _nu_plus(ombar, qbar)

    prefac = inv(qbar)

    res = -prefac * (_general_integral(num, bbar) - _general_integral(nup, bbar))

    return res
end

## Imaginary part

function _imag_ideal_dynamic_response(::NoApprox, ombar::T, qbar::T, bbar::T) where {T <: Real}

    prefac = -pi / (2 * qbar)

    num = _nu_minus(ombar, qbar)
    nup = _nu_plus(ombar, qbar)
    mubar = _chemical_potential_normalized(bbar)

    term1 = log1pexp(bbar * (mubar - num^2))
    term2 = log1pexp(bbar * (mubar - nup^2))

    return prefac * (term1 - term2) / bbar

end
