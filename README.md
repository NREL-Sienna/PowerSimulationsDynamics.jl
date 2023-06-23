# PowerSimulationsDynamics.jl

![Main - CI](https://github.com/NREL-Sienna/PowerSimulationsDynamics.jl/workflows/Main%20-%20CI/badge.svg?branch=main)
[![codecov](https://codecov.io/gh/NREL-Sienna/PowerSimulationsDynamics.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/NREL-Sienna/PowerSimulationsDynamics.jl)
[![Documentation](https://github.com/NREL-Sienna/PowerSimulationsDynamics.jl/workflows/Documentation/badge.svg)](https://nrel-sienna.github.io/PowerSimulationsDynamics.jl/stable)
[<img src="https://img.shields.io/badge/slack-@Sienna/PSID-sienna.svg?logo=slack">](https://join.slack.com/t/nrel-sienna/shared_invite/zt-glam9vdu-o8A9TwZTZqqNTKHa7q3BpQ)
[![DOI](https://zenodo.org/badge/280242020.svg)](https://zenodo.org/badge/latestdoi/280242020) [![PowerSimulationsDynamics.jl Downloads](https://shields.io/endpoint?url=https://pkgs.genieframework.com/api/v1/badge/PowerSimulationsDynamics)](https://pkgs.genieframework.com?packages=PowerSimulationsDynamics)

`PowerSimulationsDynamics.jl` is a Julia package for power system modeling and simulation of Power Systems dynamics. The objectives of the package are:

- Provide a flexible modeling framework that can accommodate different device models according to modeling needs.

- Streamline the construction of large scale differential equations problems to avoid repetition of work when adding/modifying model details.

- Exploit Julia's capabilities to improve computational performance of large scale power system dynamic simulations.

- Provide State-of-Art modeling to assess Low-Inertia Power Systems.

Check the [Project Section](https://github.com/NREL-Sienna/PowerSimulationsDynamics.jl/projects/1) to see the pipelines of new models to be added.

## Installation

```julia
julia> ]
(v1.8) pkg> add PowerSystems
(v1.8) pkg> add PowerSimulationsDynamics
```

## Usage

`PowerSimulationsDynamics.jl` uses [PowerSystems.jl](https://github.com/NREL-Sienna/PowerSystems.jl) to handle the data used in the simulations.

```julia
using PowerSimulationsDynamics
using PowerSystems
```

## Development

Contributions to the development and enahancement of PowerSimulationsDynamics is welcome. Please see [CONTRIBUTING.md](https://github.com/nrel-sienna/PowerSimulationsDynamics.jl/blob/main/CONTRIBUTING.md) for code contribution guidelines.

## License

PowerSimulationsDynamics is released under a BSD [license](https://github.com/NREL-Sienna/PowerSimulationsDynamics.jl/blob/main/LICENSE).
PowerSimulationsDynamics has been developed as part of the Scalable Integrated Infrastructure Planning (SIIP) initiative at the U.S. Department of Energy's National Renewable Energy Laboratory ([NREL](https://www.nrel.gov/))
