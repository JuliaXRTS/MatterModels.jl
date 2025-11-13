### real part

# TODO: refactor this!
function _deg_integral(x::Real)
    if abs(x) < 0.95
        # Numerically stable atanh version
        return (1 - x^2) * atanh(x)
    elseif abs(x - 1) < 1.0e-4
        # Taylor expansion around x = 1
        δ = x - 1
        return (δ^2) / 3 + (2δ^3) / 5 - (δ^4) / 7
    elseif abs(x + 1) < 1.0e-4
        # Taylor expansion around x = -1
        δ = x + 1
        return (δ^2) / 3 - (2δ^3) / 5 + (δ^4) / 7
    else
        # Use log1p if safe, otherwise fall back to log(abs(...))
        log_term = if (1 + x > 0) && (1 - x > 0)
            log1p(x) - log1p(-x)
        else
            log(abs(1 + x)) - log(abs(1 - x))
        end
        return 0.5 * (1 - x^2) * log_term
    end
end

function _real_lindhard_nonzero_temperature(::Degenerated, ombar::T, qbar::T, bbar::T) where {T <: Real}

    prefac = inv(qbar)
    num = _nu_minus(ombar, qbar)
    nup = _nu_plus(ombar, qbar)


    if qbar < 1.0e-2 * ombar
        # limiting case for large ombar
        # source: AB84, table 1, (4,2)
        prefac2 = -32.0 * qbar^2 * pi / (3 * ombar^2)
        return prefac2 * (
            1 + qbar^4 / ombar^2 + 6 * qbar^2 / (bbar * ombar^2) + 60 * qbar^4 / (bbar^2 * ombar^4)
        )
    end
    if ombar < 1.0e-4 * qbar
        # fallback on no-approx
        # TODO: find an improved formula for this case!
        return _real_lindhard_zero_temperature(NoApprox(), ombar, qbar)
    end

    # general degenerated case
    # source: AB84, eq. 10a
    return -prefac * (qbar + _deg_integral(nup) - _deg_integral(num))
end

### imaginary part

# general degenerated case
# source: AB84, eq. 26
function _imag_lindhard_nonzero_temperature(::Degenerated, ombar::T, qbar::T, bbar::T) where {T <: Real}
    prefac = -pi / (2 * qbar * bbar)

    num = _nu_minus(ombar, qbar)
    nup = _nu_plus(ombar, qbar)
    mubar = _chemical_potential_normalized(bbar)

    term1 = bbar * (one(num) - num^2)
    term2 = bbar * (one(nup) - nup^2)

    return prefac * (log1pexp(term1) - log1pexp(term2))
end
