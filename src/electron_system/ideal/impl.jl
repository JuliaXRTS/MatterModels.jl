# TODO:
# - treat ZeroTemperature mode as approximation (ZeroTemperatureApproximation)

"""
    IdealElectronSystem{T, A} <: AbstractIdealElectronSystem

Represents a non-interacting (ideal) electron system, parameterized by a numerical type (`T`)
and a response approximation model (`A`). The system may be used to model either a
finite-temperature or an effective zero-temperature electron gas. The latter is obtained
by choosing `approx = ZeroTemperatureApprox()`.

# Type Parameters
- `T <: Real`: Internal numeric type (e.g. `Float64`).
- `A <: AbstractResponseApproximation`: Approximation model used for the response
  function (e.g. `NoApprox()` or `ZeroTemperatureApprox()`).

# Fields
- `electron_density::T`: Electron number density (stored in atomic units).
- `temperature::T`: Temperature (stored in atomic units).
- `approx::A`: Approximation model for the response function.

# Constructors

    IdealElectronSystem(electron_density::T1, temperature::T2, approx::A)

Construct an ideal electron system with a given electron density, temperature,
and response approximation model.

    IdealElectronSystem(electron_density::T1, approx::ZeroTemperatureApprox)

Construct an effectively zero-temperature system.

    IdealElectronSystem(electron_density::T1, temperature::T2)

Construct an ideal electron system using the default response approximation (`NoApprox()`).

# Examples

```julia
using Unitful

# Zero-temperature system
sys1 = IdealElectronSystem(1e23u"cm^(-3)", ZeroTemperatureApprox())

# Finite-temperature system with explicit approximation
sys2 = IdealElectronSystem(1e23u"cm^(-3)", 300.0u"K", NonDegenerated())

# Finite-temperature system with default approximation
sys3 = IdealElectronSystem(1e23u"cm^(-3)", 300.0u"K")

"""
struct IdealElectronSystem{T, A <: AbstractResponseApproximation} <: AbstractIdealElectronSystem
    electron_density::T
    temperature::T
    approx::A

    function IdealElectronSystem(
            electron_density::T1,
            temperature::T2,
            approx::A
        ) where {
            T <: Real,
            T1 <: Union{T, Quantity{T}},
            T2 <: Union{T, Quantity{T}},
            A <: AbstractResponseApproximation,
        }

        ne_internal = _internalize_density(electron_density)
        temp_internal = _internalize_temperature(temperature)

        return new{T, A}(ne_internal, temp_internal)
    end
end

function IdealElectronSystem(
        electron_density::T1,
        temperature::T2
    ) where {
        T <: Real,
        T1 <: Union{T, Quantity{T}},
        T2 <: Union{T, Quantity{T}},
    }
    return IdealElectronSystem(electron_density, temperature, NoApprox())
end

function IdealElectronSystem(
        electron_density::T1,
        approx::ZeroTemperatureApprox
    ) where {
        T <: Real,
        T1 <: Union{T, Quantity{T}},
    }
    return IdealElectronSystem(electron_density, zero(T), approx)
end

@inline temperature(elsys::IdealElectronSystem) = elsys.temperature
@inline electron_density(elsys::IdealElectronSystem) = elsys.electron_density
@inline response_approximation(elsys::IdealElectronSystem) = elsys.approx
