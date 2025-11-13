# delegations

"""
    local_effective_potential(iesys::AbstractInteractingElectronSystem, om_q) -> Real

Delegates to the `local_effective_potential` method of the screening model associated with `iesys`.

# Arguments
- `iesys`: An interacting electron system.
- `om_q`: Tuple `(ω, q)` representing frequency and wavevector magnitude.

# Returns
- Local effective potential `V_eff(q, ω)` as a real number.
"""
@inline function local_effective_potential(
        iesys::AbstractInteractingElectronSystem,
        om_q::NTuple{2, T},
    ) where {T <: Real}
    return local_effective_potential(screening(iesys), om_q)
end

"""
    pseudo_potential(iesys::AbstractInteractingElectronSystem, om_q) -> Real

Delegates to the `pseudo_potential` method of the screening model associated with `iesys`.

# Arguments
- `iesys`: An interacting electron system.
- `om_q`: Tuple `(ω, q)` representing frequency and wavevector magnitude.

# Returns
- Pseudo potential `V(q, ω)` as a real number.
"""
@inline function pseudo_potential(
        iesys::AbstractInteractingElectronSystem,
        om_q::NTuple{2, T},
    ) where {T <: Real}
    return pseudo_potential(screening(iesys), om_q)
end

"""
    dielectric_function(iesys::AbstractInteractingElectronSystem, om_q) -> Real

Delegates to the `dielectric_function` using the proper electron system and screening model associated with `iesys`.

# Arguments
- `iesys`: An interacting electron system.
- `om_q`: Tuple `(ω, q)` representing frequency and wavevector magnitude.

# Returns
- Dielectric function `ϵ(q, ω)` as a real number.
"""
@inline function dielectric_function(
        iesys::AbstractInteractingElectronSystem,
        om_q::NTuple{2, T},
    ) where {T <: Real}
    return dielectric_function(proper_electron_system(iesys), screening(iesys), om_q)
end

# response functions
"""
    dynamic_response(iesys::AbstractInteractingElectronSystem, om_q) -> Complex

Compute the dynamically screened response function of the interacting electron system.

The response is computed as

```math
χ(q, ω) = \\frac{χ_0(q, ω)}{ϵ(q, ω)}
```

where χ₀ is the proper response and ϵ is the dielectric function derived from the screening model.

# Arguments

- `iesys`: An interacting electron system.
- `om_q`: Tuple (ω, q) representing frequency and wavevector magnitude.

# Returns

Complex-valued dynamic response χ(q, ω).
"""
function dynamic_response(
        iesys::AbstractInteractingElectronSystem,
        om_q::NTuple{2, T},
    ) where {T <: Real}

    # TODO: think about implementing this directly with imag and real part
    # -> this could also be used directly for imag and real part of chi_RPA
    rf = dynamic_response(proper_electron_system(iesys), om_q)
    lep = local_effective_potential(screening(iesys), om_q)
    eps = _dielectric_function(rf, lep)
    return rf / eps
end

"""
    real_dynamic_response(iesys::AbstractInteractingElectronSystem, om_q) -> Real

Return the real part of the dynamic response function of the interacting electron system.

# Returns
- Real part of `χ(q, ω)`.
"""
function real_dynamic_response(
        iesys::AbstractInteractingElectronSystem,
        om_q::NTuple{2, T},
    ) where {T <: Real}
    rf = dynamic_response(iesys, om_q)
    return real(rf)
end

"""
    imag_dynamic_response(iesys::AbstractInteractingElectronSystem, om_q) -> Real

Return the imaginary part of the dynamic response function of the interacting electron system.

# Returns
- Imaginary part of `χ(q, ω)`.
"""
function imag_dynamic_response(
        iesys::AbstractInteractingElectronSystem,
        om_q::NTuple{2, T},
    ) where {T <: Real}
    rf = dynamic_response(iesys, om_q)
    return imag(rf)
end
