@inline _fermi_wave_vector(ne::T) where {T <: Real} = cbrt(3 * pi^2 * ne)

"""
    fermi_wave_vector(esys::AbstractMatterModel) -> Real

Compute the Fermi wave vector for a given electronic system `esys`, based on its electron density.

# Arguments
- `esys`: An instance of `AbstractMatterModel` representing the electronic system.

# Returns
- Fermi wave vector `kF` as a real number.
"""
@inline fermi_wave_vector(esys::AbstractMatterModel) =
    _fermi_wave_vector(electron_density(esys))

@inline _fermi_energy_from_kF(kF::T) where {T <: Real} = kF^2 / 2

"""
    fermi_energy(esys::AbstractMatterModel) -> Real

Compute the Fermi energy of a given electronic system `esys`, using its Fermi wave vector.

# Arguments
- `esys`: An instance of `AbstractMatterModel` representing the electronic system.

# Returns
- Fermi energy `EF` as a real number.
"""
@inline fermi_energy(esys::AbstractMatterModel) =
    _fermi_energy_from_kF(fermi_wave_vector(esys))

"""
    beta(esys::AbstractMatterModel) -> Real

Compute the inverse temperature β = 1 / T of the system `esys`.

# Arguments
- `esys`: An instance of `AbstractMatterModel` representing the electronic system.

# Returns
- Inverse temperature `β` as a real number.
"""
@inline beta(esys::AbstractMatterModel) = inv(temperature(esys))

"""
    betabar(esys::AbstractMatterModel) -> Real

Compute the dimensionless inverse temperature \$\\bar\\beta= E_F / T\$ of the system `esys`.

# Arguments
- `esys`: An instance of `AbstractMatterModel` representing the electronic system.

# Returns
- Dimensionless inverse temperature `β̄` as a real number.
"""
@inline betabar(esys::AbstractMatterModel) = fermi_energy(esys) * beta(esys)

function _transform_om_q(elsys::AbstractMatterModel, om_q::NTuple{2, T}) where {T <: Real}
    kF = fermi_wave_vector(elsys)
    EF = fermi_energy(elsys)

    om, q = om_q
    ombar = om / EF
    qbar = q / kF

    return ombar, qbar
end
