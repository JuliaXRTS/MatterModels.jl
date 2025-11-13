#TBW:
# - implement an interface for ideal electron gases
# - consider adding `approximation` to the interface
# - think about generics on this interface

abstract type AbstractIdealElectronSystem <: AbstractProperElectronSystem end

"""
    response_approximation(sys::AbstractIdealElectronSystem)

Return the response approximation model associated with the given ideal electron system.
"""
function response_approximation end

"""
    _imag_ideal_dynamic_response(
        sys::AbstractIdealElectronSystem,
        approx::AbstractResponseApproximation,
        ombar::T,
        qbar::T
    ) where {T<:Real}

Compute the imaginary part of the dynamic response function of an ideal electron system
using the given response approximation. The input frequencies `ombar` and wavevectors `qbar`
are given in dimensionless units (normalized to Fermi energy and Fermi wavevector, respectively).

This is a low-level method, typically implemented by concrete approximation types.
"""
function _imag_ideal_dynamic_response end

"""

    real_ideal_dynamic_response(
        sys::AbstractIdealElectronSystem,
        approx::AbstractResponseApproximation,
        ombar::T,
        qbar::T
    ) where {T<:Real}

Compute the real part of the dynamic response function of an ideal electron system
using the given response approximation. The input frequencies `ombar` and wavevectors `qbar`
are given in dimensionless units (normalized to Fermi energy and Fermi wavevector, respectively).

This is a low-level interface meant for implementation by specific approximation models.
"""
function _real_ideal_dynamic_response end
