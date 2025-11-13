# TODO:
# - treat ZeroTemperature mode as approximation (ZeroTemperatureApproximation)

"""
    IdealElectronSystem{TM, T, A} <: AbstractIdealElectronSystem

Represents a non-interacting (ideal) electron system, parameterized by a temperature mode (`TM`),
a numerical type (`T`), and a response approximation model (`A`). This system can describe either
a zero-temperature or finite-temperature electron gas.

# Type Parameters
- `TM <: AbstractTemperatureMode`: Temperature mode (`ZeroTemperature` or `FiniteTemperature`).
- `T <: Real`: Internal numeric type (e.g., `Float64`).
- `A <: AbstractResponseApproximation`: Approximation used for the response function.

# Fields
- `electron_density::T`: Electron number density (stored in atomic units).
- `temperature::T`: Temperature (stored in atomic units).
- `approx::A`: Approximation model for the response function (e.g., `NoApprox()`).

# Temperature Modes

- **ZeroTemperature**: Models a zero-temperature Fermi gas. Temperature is internally set to zero.
- **FiniteTemperature**: Models finite-temperature effects in the electron gas.

# Constructors

    IdealElectronSystem{TM}(electron_density::T1, temperature::T2, approx::A)

Construct an ideal electron system with a specified temperature mode, electron density, temperature,
and approximation model.

    IdealElectronSystem{ZeroTemperature}(electron_density::T1, approx::A)

Construct a zero-temperature system with a custom approximation.

    IdealElectronSystem{ZeroTemperature}(electron_density::T1)

Construct a zero-temperature system using the default approximation (`NoApprox()`).

    IdealElectronSystem(electron_density::T1, temperature::T2, approx::A)

Construct a finite-temperature system. The temperature mode is inferred as `FiniteTemperature`.

    IdealElectronSystem{FiniteTemperature}(electron_density::T1, temperature::T2)

Construct a finite-temperature system with the default approximation (`NoApprox()`).

    IdealElectronSystem(electron_density::T1, temperature::T2)

Alias for `IdealElectronSystem{FiniteTemperature}(...)`.

# Examples

```julia
using Unitful

# Zero-temperature system with default approximation
sys1 = IdealElectronSystem{ZeroTemperature}(1e23u"cm^(-3)")

# Finite-temperature system with explicit approximation
sys2 = IdealElectronSystem(1e23u"cm^(-3)", 300.0u"K", NonDegenerated())

# Finite-temperature system with default approximation
sys3 = IdealElectronSystem(1e23u"cm^(-3)", 300.0u"K")
```
"""
struct IdealElectronSystem{TM, T, A <: AbstractResponseApproximation} <: AbstractIdealElectronSystem
    electron_density::T
    temperature::T
    approx::A

    function IdealElectronSystem{TM}(
            electron_density::T1,
            temperature::T2,
            approx::A
        ) where {
            T <: Real,
            T1 <: Union{T, Quantity{T}},
            T2 <: Union{T, Quantity{T}},
            TM <: AbstractTemperatureMode,
            A <: AbstractResponseApproximation,
        }

        ne_internal = _internalize_density(electron_density)
        temp_internal = _internalize_temperature(temperature)

        return new{TM, T, A}(ne_internal, temp_internal)
    end
end

# TODO: think about making this a default constructor
function IdealElectronSystem{ZeroTemperature}(
        electron_density::T1,
        approx::A
    ) where {
        T <: Real,
        T1 <: Union{T, Quantity{T}},
        A <: AbstractZeroTemperatureApproximation,
    }
    return IdealElectronSystem{ZeroTemperature}(electron_density, zero(T), approx)
end

function IdealElectronSystem{ZeroTemperature}(
        electron_density::T1,
    ) where {
        T <: Real,
        T1 <: Union{T, Quantity{T}},
    }
    return IdealElectronSystem{ZeroTemperature}(electron_density, zero(T), NoApprox())
end

function IdealElectronSystem(
        electron_density::T1,
        temperature::T2,
        approx::A
    ) where {
        T <: Real,
        T1 <: Union{T, Quantity{T}},
        T2 <: Union{T, Quantity{T}},
        A <: AbstractFiniteTemperatureApproximation,
    }
    return IdealElectronSystem{FiniteTemperature}(electron_density, temperature, approx)
end

function IdealElectronSystem{FiniteTemperature}(
        electron_density::T1,
        temperature::T2,
    ) where {
        T <: Real,
        T1 <: Union{T, Quantity{T}},
        T2 <: Union{T, Quantity{T}},
    }
    return IdealElectronSystem{FiniteTemperature}(electron_density, temperature, NoApprox())
end

function IdealElectronSystem(
        electron_density::T1,
        temperature::T2,
    ) where {
        T <: Real,
        T1 <: Union{T, Quantity{T}},
        T2 <: Union{T, Quantity{T}},
    }
    return IdealElectronSystem{FiniteTemperature}(electron_density, temperature, NoApprox())
end

@inline temperature(elsys::IdealElectronSystem{ZeroTemperature, T}) where {T <: Real} = zero(T)
@inline temperature(elsys::IdealElectronSystem{FiniteTemperature, T}) where {T <: Real} =
    elsys.temperature
@inline electron_density(elsys::IdealElectronSystem) = elsys.electron_density
@inline response_approximation(elsys::IdealElectronSystem) = elsys.approx

### default implementations -> needed for separation between zero and nonzero case

# Zero temperature case

function _imag_ideal_dynamic_response(
        elsys::IdealElectronSystem{ZeroTemperature},
        approx::AbstractZeroTemperatureApproximation,
        ombar::T,
        qbar::T
    ) where {T <: Real}


    return _imag_lindhard_zero_temperature(approx, ombar, qbar)
end

function _real_ideal_dynamic_response(
        elsys::IdealElectronSystem{ZeroTemperature},
        approx::AbstractZeroTemperatureApproximation,
        ombar::T,
        qbar::T
    ) where {T <: Real}

    return _real_lindhard_zero_temperature(approx, ombar, qbar)
end


# finite temperature case

function _imag_ideal_dynamic_response(
        elsys::IdealElectronSystem{FiniteTemperature},
        approx::AbstractFiniteTemperatureApproximation,
        ombar::T,
        qbar::T
    ) where {T <: Real}

    return _imag_lindhard_nonzero_temperature(approx, ombar, qbar, betabar(elsys))
end

function _real_ideal_dynamic_response(
        elsys::IdealElectronSystem{FiniteTemperature},
        approx::AbstractFiniteTemperatureApproximation,
        ombar::T,
        qbar::T
    ) where {T <: Real}
    return _real_lindhard_nonzero_temperature(approx, ombar, qbar, betabar(elsys))
end
