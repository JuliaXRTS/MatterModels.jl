function _nu_minus(ombar, qbar)
    return (ombar - qbar^2) / (2 * qbar)
end

function _nu_plus(ombar, qbar)
    return (ombar + qbar^2) / (2 * qbar)
end

function heaviside(t)
    return 0.5 * (sign(t) + 1)
end

# Ichimaru:2018
function _chemical_potential_normalized(bbar)
    A = 0.25954
    B = 0.072
    b = 0.858

    term1 = 3 / 2 * log(bbar)
    term2 = log(4 / (3 * sqrt(pi)))
    term3 = (A * bbar^(b + 1) + B * bbar^((b + 1) / 2)) / (1 + A * bbar^b)

    return (term1 + term2 + term3) / bbar
end

# GV:2005 eq. (4.43)
function _F(x, bbar)
    mubar = _chemical_potential_normalized(bbar)
    denom = exp(bbar * (x^2 - mubar)) + one(bbar)
    return x / denom
end


"""

    _stable_log_term(x::Real)

Stable version of

```math
\\log\\left(\\left\\vert\\frac{1+x}{1-x}\\right\\vert\\right)
```

"""
function _stable_log_term(x::Real)
    return abs(x) < 1 ? 2 * atanh(x) : 2 * atanh(inv(x))
end
