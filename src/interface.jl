"""
    AbstractMatterModel

Abstract base type representing a generic matter model.

Concrete subtypes should implement interface functions such as `temperature` and `electron_density`
to provide physical properties of the matter system.
"""
abstract type AbstractMatterModel end

"""
    temperature(m::AbstractMatterModel)

Return the temperature of the given matter model `m` in internal units.

This is an interface function that must be implemented by all concrete matter models.
"""
function temperature end

"""
    electron_density(m::AbstractMatterModel)

Return the electron density of the given matter model `m`.

This is an interface function that must be implemented by all concrete matter models.
"""
function electron_density end

# TODO:
# - consider moving `electron_density` downstream, maybe to the interface for electron
# systems
# - implement `temperature(::Unit,::AbstactMatterModel)` to return the temperature in a
# user defined unit; similar for `electron_density`
# - propagate the unitful call downstream, e.g. to `beta`
