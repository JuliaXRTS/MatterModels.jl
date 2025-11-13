# local field correction
"""
    NoLocalFieldCorrection <: AbstractLocalFieldCorrection

Represents a model with no local field correction (i.e., the correction is zero).
"""
struct NoLocalFieldCorrection <: AbstractLocalFieldCorrection end

function _compute(::NoLocalFieldCorrection, om_q::NTuple{2, T}) where {T <: Real}
    return zero(T)
end

# pseudo potential
"""
    CoulombPseudoPotential <: AbstractPseudoPotential

Represents the bare Coulomb pseudo potential: V(q) = e² / q².
"""
struct CoulombPseudoPotential <: AbstractPseudoPotential end

# TODO: insert 4pi!
function _compute(::CoulombPseudoPotential, om_q::NTuple{2, T}) where {T <: Real}
    om, q = om_q
    return ELEMENTARY_CHARGE_SQUARED / q^2
end

### screening
# no screening
"""
    NoScreening <: AbstractScreening

A screening model with no dielectric screening; dielectric function is unity, and corrections vanish.
"""
struct NoScreening <: AbstractScreening end

"""
    dielectric_function(esys, scr::NoScreening, om_q) -> Real

Return the dielectric function for a `NoScreening` model, which is always 1.
"""
dielectric_function(
    esys::AbstractProperElectronSystem,
    scr::NoScreening,
    om_q::NTuple{2, T},
) where {T <: Real} = one(T)

"""
    local_effective_potential(scr::NoScreening, om_q) -> Real

Return the local effective potential for the `NoScreening` model, which is zero.
"""
local_effective_potential(src::NoScreening, om_q::NTuple{2, T}) where {T <: Real} = zero(T)

"""
    pseudo_potential(scr::NoScreening, om_q) -> Real

Return the pseudo potential for the `NoScreening` model, which is zero.
"""
@inline pseudo_potential(scr::NoScreening, om_q::NTuple{2, T}) where {T <: Real} = zero(T)

"""
    local_field_correction(scr::NoScreening, om_q) -> Real

Return the local field correction for the `NoScreening` model, which is zero.
"""
@inline local_field_correction(scr::NoScreening, om_q::NTuple{2, T}) where {T <: Real} =
    zero(T)

# screeing with pseudo potential and local field correction
"""
    Screening{PP, LFC} <: AbstractScreening

A general screening model that combines a pseudo potential and a local field correction.

# Fields
- `pseudo_potential`: An instance of `AbstractPseudoPotential`.
- `lfc`: An instance of `AbstractLocalFieldCorrection`.
"""
struct Screening{PP <: AbstractPseudoPotential, LFC <: AbstractLocalFieldCorrection} <:
    AbstractScreening
    pseudo_potential::PP
    lfc::LFC
end

Screening() = Screening(CoulombPseudoPotential(), NoLocalFieldCorrection())

@inline pseudo_potential(scr::Screening, om_q) = _compute(scr.pseudo_potential, om_q)
@inline local_field_correction(scr::Screening, om_q) = _compute(scr.lfc, om_q)
