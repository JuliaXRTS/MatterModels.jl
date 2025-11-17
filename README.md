# ElectronicStructureModels.jl

[![Stable Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaXRTS.github.io/ElectronicStructureModels.jl/stable)
[![Development documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuliaXRTS.github.io/ElectronicStructureModels.jl/dev)
[![Test workflow status](https://github.com/JuliaXRTS/ElectronicStructureModels.jl/actions/workflows/Test.yml/badge.svg?branch=main)](https://github.com/JuliaXRTS/ElectronicStructureModels.jl/actions/workflows/Test.yml?query=branch%3Amain)
[![BestieTemplate](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/JuliaBesties/BestieTemplate.jl/main/docs/src/assets/badge.json)](https://github.com/JuliaBesties/BestieTemplate.jl)
[![code style: runic](https://img.shields.io/badge/code_style-%E1%9A%B1%E1%9A%A2%E1%9A%BE%E1%9B%81%E1%9A%B2-black)](https://github.com/fredrikekre/Runic.jl)

## Overview

**ElectronicStructureModels.jl** is an open-source Julia package providing models from
electronic-structure theory for use in larger simulations and data analysis workflows. It
is designed with a focus on applications in **X-ray Thomson Scattering (XRTS)**, but the
components are general enough to be used in a broader context involving electron-gas
models and response functions.

## Key Features

The package currently includes:

- Dynamic structure factors (DSF) for the uniform electron gas; including the finite-temperature Lindhard function within the Random Phase Approximation (RPA)
- Several approximation levels for evaluating DSF and related quantities
- A clean, extensible interface for adding new electronic-structure models in the future

## Installation

The package is not registered yet. To install it locally:

```bash
git clone https://github.com/JuliaXRTS/ElectronicStructureModels.jl.git
cd ElectronicStructureModels.jl
julia --project -e "import Pkg; Pkg.build()"
```

Then in Julia:

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
```

This will also install required dependencies such as `Unitful.jl`.

## Quickstart Example

Below is a minimal example showing how to construct an electron system and compute a
dynamic structure factor. The example uses `CairoMakie` for plotting, but any plotting
backend can be used.

```julia
using ElectronicStructureModels
using Unitful
using CairoMakie

NE   = 10e23u"cm^(-3)"   # electron density
TEMP = 1.0u"eV"          # temperature

ideal = IdealElectronSystem(NE, TEMP) # idear proper system
esys  = InteractingElectronSystem(ideal, Screening()) # RPA screening model

KF = fermi_wave_vector(esys) # Fermi wave vector (for scaling)
EF = fermi_energy(esys) # Fermi energy (for scaling)

ω_arr = range(1e-4, 6, 100) .* EF # array of energy transfer values
q     = 1.6 * KF # momentum transfer value

# calculation of the DSF
dsf_values = [dynamic_structure_factor(esys, (ω, q)) for ω in ω_arr]


# Plotting of the DSF as a function of the energy transfer normalized to unit height.
fig = Figure()
ax  = Axis(fig[1,1])

norm = maximum(dsf_values)
lines!(ax, ω_arr ./ EF, dsf_values ./ norm)

display(fig)
```

## Contributing

Contributions of all kinds are welcome.
Please read the [contributing guide](docs/src/90-contributing.md) or its rendered version in the [online documentation](https://JuliaXRTS.github.io/ElectronicStructureModels.jl/dev/90-contributing/) before opening pull requests or issues.

## Credits and Funding

This work has been partially supported by the **Center for Advanced Systems Understanding (CASUS)**, funded by the German Federal Ministry of Education and Research (BMBF) and the Saxon Ministry for Science, Culture and Tourism (SMWK) with tax funds based on the budget approved by the Saxon State Parliament.

### Acknowledgements

Special thanks to **Tobias Dornheim** for funding, valuable input and support.

## License

[MIT](LICENSE) © Uwe Hernandez Acosta
