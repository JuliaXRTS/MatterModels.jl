### general electron systems


"""
    AbstractElectronSystem <: AbstractMatterModel

Abstract base type for electron matter models.

## Interface requirements

Concrete subtypes must implement the following functions:

1. Matter model interface:
- `electron_density(sys::AbstractElectronSystem)::Real`
- `temperature(sys::AbstractElectronSystem)::Real`

2. Electronic system interface:
- `imag_dynamic_response(sys::AbstractElectronSystem, om_q::NTuple{2, Real})::Real`
- `real_dynamic_response(sys::AbstractElectronSystem, om_q::NTuple{2, Real})::Real`
"""
abstract type AbstractElectronSystem <: AbstractMatterModel end

"""
    imag_dynamic_response(sys::AbstractElectronSystem, om_q::NTuple{2, T}) where {T<:Real}

Return the imaginary part of the dynamic response function for the given electron system `sys`
at frequency and wavevector specified by `om_q = (ω, q)`.
"""
function imag_dynamic_response end

"""
    real_dynamic_response(sys::AbstractElectronSystem, om_q::NTuple{2, T}) where {T<:Real}

Return the real part of the dynamic response function for the given electron system `sys`
at frequency and wavevector specified by `om_q = (ω, q)`.
"""
function real_dynamic_response end

### proper (non-interacting) electron system
# - electron system without screening

abstract type AbstractProperElectronSystem <: AbstractElectronSystem end
