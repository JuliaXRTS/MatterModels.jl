module MatterModels

hello_world() = "Hello, World!"

# temperature
export FiniteTemperature, ZeroTemperature

# matter model
export AbstractMatterModel

# Electron system
export AbstractElectronSystem
export temperature, electron_density, imag_dynamic_response, real_dynamic_response
export fermi_wave_vector,
    fermi_energy, beta, betabar, dynamic_response, dynamic_structure_factor

export AbstractProperElectronSystem
export AbstractInteractingElectronSystem
export proper_electron_system, screening

# screening
export AbstractScreening, Screening, NoScreening
export dielectric_function,
    pseudo_potential, local_field_correction, local_effective_potential
export AbstractPseudoPotential, CoulombPseudoPotential
export AbstractLocalFieldCorrection, NoLocalFieldCorrection

# concrete electron systems
export IdealElectronSystem
export AbstractResponseApproximation, NoApprox, NonDegenerated, Degenerated
export response_approximation
export InteractingElectronSystem

# constants
export HBARC,
    HBARC_eV_ANG,
    ELECTRONMASS,
    ALPHA,
    ALPHA_SQUARE,
    ELEMENTARY_CHARGE_SQUARED,
    ELEMENTARY_CHARGE,
    HARTREE,
    BOHR_RADIUS_ANG

using QEDcore
using QuadGK
using Unitful
using LogExpFunctions
using SpecialFunctions

include("utils.jl")
include("units.jl")
include("constants.jl")
include("temperature.jl")
include("interface.jl")
include("generic.jl")
include("lookup.jl")
include("data_driven/impl.jl")
include("electron_system/utils.jl")
include("electron_system/interface.jl")
include("electron_system/generic.jl")
include("electron_system/ideal/approximations/interface.jl")
include("electron_system/ideal/approximations/no_approx.jl")
include("electron_system/ideal/approximations/non_degenerated.jl")
include("electron_system/ideal/approximations/degenerated.jl")
include("electron_system/ideal/utils.jl")
include("electron_system/ideal/interface.jl")
include("electron_system/ideal/generic.jl")
include("electron_system/ideal/impl.jl")
include("electron_system/interacting/screening/interface.jl")
include("electron_system/interacting/screening/generic.jl")
include("electron_system/interacting/screening/impl.jl")
include("electron_system/interacting/interface.jl")
include("electron_system/interacting/generic.jl")
include("electron_system/interacting/impl.jl")

end
