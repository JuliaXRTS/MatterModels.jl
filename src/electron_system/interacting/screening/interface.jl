"""

    AbstractScreening

Abstract base type for screening models. Interface functions, which must be implemented:

- `pseudo_potential(::AbstractInteractingElectronSystem)::AbstractPseudoPotential`
- `local_field_correction(::AbstractInteractingElectronSystem)::AbstractLocalFieldCorrection`

"""
abstract type AbstractScreening end
abstract type AbstractPseudoPotential end
abstract type AbstractLocalFieldCorrection end

"""

    pseudo_potential(::AbstractScreening, om_q::Tuple)

Interface function: return value of the pseudo potential for interacting electron systems
"""
function pseudo_potential end

"""

    local_field_correction(::AbstractScreening, om_q::Tuple)

Interface function: return value of the local field correction for interacting electron systems
"""
function local_field_correction end
