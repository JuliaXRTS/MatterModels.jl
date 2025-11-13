# local effective potential
@inline function _local_effective_potential(vp, lfc)
    return vp * (one(lfc) - lfc)
end

function local_effective_potential(
        scr::AbstractScreening,
        om_q::NTuple{2, T},
    ) where {T <: Real}
    vp = pseudo_potential(scr, om_q)
    lfc = local_field_correction(scr, om_q)

    return _local_effective_potential(vp, lfc)
end

"""

    dielectric_function(::AbstractIdealElectronSystem,::AbstractScreening,om_q)

Return the value of the dielectric function for the given proper response function and
given local effective potential.
"""
function _dielectric_function(rf::Number, lep::Number)
    temp = lep * rf
    return one(temp) - temp
end

"""

    dielectric_function(::AbstractIdealElectronSystem,::AbstractScreening,om_q)

Return the value of the dielectric function for the given electron system and given screening.
"""
function dielectric_function(
        esys::AbstractProperElectronSystem,
        scr::AbstractScreening,
        om_q::NTuple{2, T},
    ) where {T <: Real}
    rf = dynamic_response(esys, om_q)
    lep = local_effective_potential(scr, om_q)
    return _dielectric_function(rf, lep)
end
