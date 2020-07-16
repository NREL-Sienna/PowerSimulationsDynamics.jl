# PowerSimulations

[![Build Status](https://img.shields.io/travis/com/NREL-SIIP/PowerSimulationsDynamic.jl/master.svg)](https://travis-ci.com/NREL-SIIP/PowerSimulationsDynamic.jl)
[![codecov](https://codecov.io/gh/NREL-SIIP/PowerSimulationsDynamic.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/NREL-SIIP/PowerSimulationsDynamic.jl)
[![Documentation](https://github.com/NREL-SIIP/PowerSimulationsDynamic.jl/workflows/Documentation/badge.svg)](https://nrel-siip.github.io/PowerSimulationsDynamic.jl/latest)
[![Join the chat at https://gitter.im/NREL/PowerSimulationsDynamic.jl](https://badges.gitter.im/NREL/PowerSimulationsDynamic.jl.svg)](https://gitter.im/NREL/PowerSimulationsDynamic.jl)

## The current implementation of the functionalities can be seen in the test codes.

`PowerSimulationsDynamic.jl` is a Julia package for power system modeling and simulation of Power Systems dynamics. The objectives of the package are:

- Provide a flexible modeling framework that can accommodate different device models according to modeling needs.

- Streamline the construction of large scale differential equations problems to avoid repetition of work when adding/modifying model details.

- Exploit Julia's capabilities to improve computational performance of large scale power system dynamic simulations.

- Provide State-of-Art modeling to assess Low-Innertia Power Systems.

## Installation

```julia
julia> ]
(v1.3) pkg> add PowerSystems
(v1.3) pkg> add PowerSimulationsDynamic
```
## Usage

`PowerSimulationsDynamic.jl` uses [PowerSystems.jl](https://github.com/NREL/PowerSystems.jl) to handle the data used in the simulations.

```julia
using PowerSimulations
using PowerSystems
```

## Development

Contributions to the development and enahancement of PowerSimulationsDynamic is welcome. Please see [CONTRIBUTING.md](https://github.com/NREL-SIIP/PowerSimulationsDynamic.jl/blob/master/CONTRIBUTING.md) for code contribution guidelines.

## License

PowerSimulationsDynamic is released under a BSD [license](https://github.com/NREL/PowerSimulationsDynamic.jl/blob/master/LICENSE). PowerSimulations has been developed as part of the Scalable Integrated Infrastructure Planning (SIIP)
initiative at the U.S. Department of Energy's National Renewable Energy Laboratory ([NREL](https://www.nrel.gov/))
