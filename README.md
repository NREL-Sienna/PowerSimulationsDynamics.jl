# PowerSimulationsDynamics.jl

[![Build Status](https://travis-ci.com/NREL-SIIP/PowerSimulationsDynamics.jl.svg?branch=master)](https://travis-ci.com/NREL-SIIP/PowerSimulationsDynamics.jl)
[![codecov](https://codecov.io/gh/NREL-SIIP/PowerSimulationsDynamics.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/NREL-SIIP/PowerSimulationsDynamics.jl)
[![Documentation](https://github.com/NREL-SIIP/PowerSimulationsDynamics.jl/workflows/Documentation/badge.svg)](https://nrel-siip.github.io/PowerSimulationsDynamics.jl/latest)
[![Join the chat at https://gitter.im/NREL/PowerSimulationsDynamics.jl](https://badges.gitter.im/NREL/PowerSimulationsDynamics.jl.svg)](https://gitter.im/NREL/PowerSimulationsDynamics.jl)

`PowerSimulationsDynamics.jl` is a Julia package for power system modeling and simulation of Power Systems dynamics. The objectives of the package are:

- Provide a flexible modeling framework that can accommodate different device models according to modeling needs.

- Streamline the construction of large scale differential equations problems to avoid repetition of work when adding/modifying model details.

- Exploit Julia's capabilities to improve computational performance of large scale power system dynamic simulations.

- Provide State-of-Art modeling to assess Low-Innertia Power Systems.

## Installation

```julia
julia> ]
(v1.3) pkg> add PowerSystems
(v1.3) pkg> add PowerSimulationsDynamics
```
## Usage

`PowerSimulationsDynamics.jl` uses [PowerSystems.jl](https://github.com/NREL/PowerSystems.jl) to handle the data used in the simulations.

```julia
using PowerSimulationsDynamics
using PowerSystems
```

## Development

Contributions to the development and enahancement of PowerSimulationsDynamic is welcome. Please see [CONTRIBUTING.md](https://github.com/NREL-SIIP/PowerSimulationsDynamics.jl/blob/master/CONTRIBUTING.md) for code contribution guidelines.

## License

PowerSimulationsDynamics is released under a BSD [license](https://github.com/NREL/PowerSimulationsDynamics.jl/blob/master/LICENSE). PowerSimulations has been developed as part of the Scalable Integrated Infrastructure Planning (SIIP)
initiative at the U.S. Department of Energy's National Renewable Energy Laboratory ([NREL](https://www.nrel.gov/))
