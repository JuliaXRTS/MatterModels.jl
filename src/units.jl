# TODO:
# - refactor this!
# - add more units with energy dimention (current, time, momentum?)
# - move this to `NaturallyUnitful` (`naturally` is type unstable!)

const speed_of_light = 1u"c"
const hbar = Unitful.ħ
const kBoltzmann = Unitful.k
const hbarc = hbar * speed_of_light

function _get_dim_power(x::Unitful.Dimensions{T}) where {T}
    return prod(_get_dim_power.(T))
end

_get_dim_power(x::Unitful.Dimension{:Length}) = -x.power
_get_dim_power(x::Unitful.Dimension{:Time}) = -x.power
_get_dim_power(x::Unitful.Dimension{:Mass}) = x.power
_get_dim_power(x::Unitful.Dimension{:Temperature}) = x.power
#_get_dim_power(x::Unitful.Dimension{:Current}) = x.power

function _get_dim_factor(x::Unitful.Dimensions{T}) where {T}
    return prod(@. _get_dim_factor(T, _get_dim_power(T)))
end

function _get_dim_factor(x::Unitful.Dimension{:Length}, p::Rational)
    return hbarc^p
end

function _get_dim_factor(x::Unitful.Dimension{:Mass}, p::Rational)
    return speed_of_light^(2 * p)
end

function _get_dim_factor(x::Unitful.Dimension{:Temperature}, p::Rational)
    return kBoltzmann^p
end


function mynatural(q; base = u"eV")
    d = dimension(q)
    p = _get_dim_power(d)
    nunit = base^(p)
    nfac = _get_dim_factor(d)
    return uconvert(nunit, q * nfac)
end

# opt out, if the quanity is already in eV
gettypeparams(::Unitful.FreeUnits{T, U, V}) where {T, U, V} = T, U, V
const eV_DIM = gettypeparams(u"GeV")[2]

function mynatural(q::Quantity{T, eV_DIM}; base = u"eV") where {T <: Number}
    return uconvert(base, q)
end

# FIXME: generalize to powers `eV^p`, currently it errors for `p!=1

## conversion for structure factors

const me = ustrip(mynatural(Unitful.me))

function _eV2me_dimless(t_eV, p::Int = 1)
    return ustrip(t_eV) * me^p
end

function _density2me_dimless(density)
    density_eV = mynatural(density)
    return _eV2me_dimless(density_eV, -3)
end

function _energy2me_dimless(energy)
    energy_eV = mynatural(energy)
    return _eV2me_dimless(energy_eV, -1)
end

# simple unit transforms
# TODO:
# - write unit tests for this
# - if other unit transforms pop up, add them here
#
_bohr2ang(x_bohr, p = 1) = x_bohr * (BOHR_RADIUS_ANG^p)
_bohr2cm(x_bohr, p = 1) = x_bohr * ((BOHR_RADIUS_ANG * 1.0e-8)^p)
_bohr2fm(x_bohr, p = 1) = x_bohr * ((BOHR_RADIUS_ANG * 1.0e5)^p)
_bohr2inv_eV(x_bohr, p = 1) = x_bohr * ((BOHR_RADIUS_ANG / HBARC_eV_ANG)^p)
_bohr2inv_me(x_bohr, p = 1) = x_bohr * ((BOHR_RADIUS_ANG / HBARC_eV_ANG * ELECTRONMASS)^p)

_inv_bohr2ang(x_inv_bohr, p = 1) = x_inv_bohr / (BOHR_RADIUS_ANG^p)
_inv_bohr2cm(x_inv_bohr, p = 1) = x_inv_bohr / ((BOHR_RADIUS_ANG * 1.0e-8)^p)
_inv_bohr2fm(x_inv_bohr, p = 1) = x_inv_bohr / ((BOHR_RADIUS_ANG * 1.0e5)^p)
_inv_bohr2eV(x_inv_bohr, p = 1) = x_inv_bohr / ((BOHR_RADIUS_ANG / HBARC_eV_ANG)^p)
_inv_bohr2me(x_inv_bohr, p = 1) = x_inv_bohr / ((BOHR_RADIUS_ANG / HBARC_eV_ANG * ELECTRONMASS)^p)

_ang2inv_eV(x_ang, p = 1) = x_ang / (HBARC_eV_ANG^p)
_inv_ang2eV(x_inv_ang, p = 1) = x_inv_ang * (HBARC_eV_ANG^p)

_me2eV(x_me, p = 1) = x_me * (ELECTRONMASS)^p
_me2keV(x_me, p = 1) = x_me * (ELECTRONMASS * 1.0e-3)^p
_me2MeV(x_me, p = 1) = x_me * (ELECTRONMASS * 1.0e-6)^p
_me2hartree(x_me, p = 1) = x_me * (ELECTRONMASS / HARTREE)^p
_me2inv_ang(x_me, p = 1) = x_me * (ELECTRONMASS / HBARC_eV_ANG)^p
_me2inv_bohr(x_me, p = 1) = x_me * (ELECTRONMASS / HBARC_eV_ANG * BOHR_RADIUS_ANG)^p

# this is inv(Ang) to eV!
function _convert_Ang_to_eV(q)
    return q * HBARC_eV_ANG
end

# this is eV to inv Ang!
function _convert_eV_to_Ang(q)
    return q / HBARC_eV_ANG
end
