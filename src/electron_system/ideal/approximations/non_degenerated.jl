function _non_deg_integal(nu, beta)
    mu = _chemical_potential_normalized(beta)
    prefac = -exp(mu * beta) / beta

    return prefac * dawson(nu * sqrt(beta)) * sqrt(pi)
end

function _real_lindhard_nonzero_temperature(::NonDegenerated, ombar::T, qbar::T, bbar::T) where {T <: Real}
    num = _nu_minus(ombar, qbar)
    nup = _nu_plus(ombar, qbar)

    prefac = inv(qbar)
    return -prefac * (_non_deg_integal(num, bbar) - _non_deg_integal(nup, bbar))
end


function _imag_lindhard_nonzero_temperature(::NonDegenerated, ombar::T, qbar::T, bbar::T) where {T <: Real}
    prefac = -pi / (2 * qbar * bbar)

    num = _nu_minus(ombar, qbar)
    nup = _nu_plus(ombar, qbar)
    mubar = _chemical_potential_normalized(bbar)

    term1 = bbar * (mubar - num^2)
    term2 = bbar * (mubar - nup^2)

    return prefac * (exp(term1) - exp(term2))
end
