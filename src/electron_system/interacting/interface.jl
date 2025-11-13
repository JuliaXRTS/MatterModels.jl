"""

    AbstractInteractingElectronSystem

Abstract base type for electronic system with interacting electrons by introducing screening potentials.
Interface functions to be implemented are:

1. Matter model interface:
    - `electron_density(::AbstractElectronSystem)::Real`
    - `temperature(::AbstractElectronSystem)::Real`

2. Electronic system interface:
    - `imag_dynamic_response(::AbstractElectronSystem, om_q::NTuple{2,Real})::Real`
    - `real_dynamic_response(::AbstractElectronSystem, om_q::

3. Interacting electronic system interface:
    - `proper_electron_system(::AbstractInteractingElectronSystem)::AbstractProperElectronSystem`
    - `screening(::AbstractInteractingElectronSystem)::AbstractScreening`
"""
abstract type AbstractInteractingElectronSystem <: AbstractElectronSystem end

"""

    proper_electron_system(::AbstractInteractingElectronSystem)

Interface function: return the proper electron system associated with the interacting system.
"""
function proper_electron_system end

"""

    screening(::AbstractInteractingElectronSystem)

Interface function: return the screening associated with the interacting system.
"""
function screening end
