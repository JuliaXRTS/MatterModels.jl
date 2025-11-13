# Temperature modes
abstract type AbstractTemperatureMode end
struct FiniteTemperature <: AbstractTemperatureMode end
struct ZeroTemperature <: AbstractTemperatureMode end
