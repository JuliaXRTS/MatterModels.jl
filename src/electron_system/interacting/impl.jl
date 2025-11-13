"""
    InteractingElectronSystem{PES, S} <: AbstractInteractingElectronSystem

Represents an interacting electron system consisting of a proper electron system (e.g., a free or ideal Fermi gas)
and a screening model that modifies its response.

# Type Parameters
- `PES <: AbstractProperElectronSystem`: Type of the proper (non-interacting) electron system.
- `S <: AbstractScreening`: Type of the screening model.

# Fields
- `proper_elec_system::PES`: The underlying proper electron system.
- `screening::S`: The screening model used to compute interactions.

# Constructors

```julia
InteractingElectronSystem(pelsys::PES, scr::S)
```

Construct an interacting electron system directly from a proper electron system and a screening model.

```julia
InteractingElectronSystem{PES}(electron_density::T1, temp::T2, scr::S)
```
Construct an interacting system by providing physical parameters, where `PES` is the desired proper electron system type.

# Examples

```julia
iesys = InteractingElectronSystem{IdealFermiGas}(1e23, 300.0, Screening())
```
"""
struct InteractingElectronSystem{PES <: AbstractProperElectronSystem, S <: AbstractScreening} <:
    AbstractInteractingElectronSystem
    proper_elec_system::PES
    screening::S

    function InteractingElectronSystem(
            pelsys::PES,
            scr::S,
        ) where {PES <: AbstractProperElectronSystem, S <: AbstractScreening}

        return new{PES, S}(pelsys, scr)
    end

    function InteractingElectronSystem{PES}(
            electron_density::T1,
            temp::T2,
            scr::S,
        ) where {
            T <: Real,
            T1 <: Union{T, Quantity},
            T2 <: Union{T, Quantity},
            PES <: AbstractProperElectronSystem,
            S <: AbstractScreening,
        }
        proper_elec_system = PES(electron_density, temp)
        return new{PES, S}(proper_elec_system, scr)
    end
end

# accessors
@inline proper_electron_system(elsys::InteractingElectronSystem) = elsys.proper_elec_system
@inline screening(elsys::InteractingElectronSystem) = elsys.screening
@inline temperature(elsys::InteractingElectronSystem) =
    temperature(proper_electron_system(elsys))
@inline electron_density(elsys::InteractingElectronSystem) =
    electron_density(proper_electron_system(elsys))
