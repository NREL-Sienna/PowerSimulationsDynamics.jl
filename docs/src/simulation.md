# Simulation Structure

After constructing the System data from `PowerSystems.jl` with its dynamic components, a `Simulation` structure must be constructed. To do so, the arguments given are the following:

```raw
Simulation(
    Type of Model,
    System,
    Folder Directory,
    Time Span,
    Optional: Perturbation,
    Optional: Keyword Arguments,
)
```

The Simulation will construct all the indexing for the states and will attempt to initialize the dynamical system in steady state based on the power flow data provided.

We will go through each argument that must be included:

## Type of Simulation Model

In `PowerSimulationDynamics` there are two types of models depending on which Solver will be used to execute the Simulation. More details are discussed in [Models](models.md) subsection:

- The **ResidualModel** is used when a residual model of the type ``res = f(x) - M\cdot dx`` is constructed. This allows to use pure DAE solvers, such as Sundials IDA.
- The **MassMatrixModel** is used when a mass matrix model of the type ``M \cdot dx = f(x)`` is constructed. This allows to use mass matrix ODE solvers, such as Rodas4 or Rodas5.

## System

The data system constructed using `PowerSystems.jl`. At least one dynamic device must be included in the system to construct a Simulation.

## Folder Directory

The directory that can be used to serialize the system if required. Usually `pwd()`, the current directory, is used.

## Time Span

A tuple of `Float` must be used to define the time span (in seconds) for the Simulation. For example `(0.0, 20.0)` is a valid simulation that starts at ``t = 0`` seconds that finishes at ``t = 20`` seconds.

## Optional: Perturbations

A perturbation structure must be defined that will alter the system from its steady state operation. More details on the available perturbations in their respective [Section](perturbations.md).

## Optional: Keyword Arguments

Different keyword arguments can be included in the Simulation for different purposes


### `initialize_simulation::Bool`

Default: `true`. Runs the initialization routine. If set to `false`, it will use the operating point for voltages stored in the System. Dynamic states will be initialized at zero.

### `initial_condition::Vector{Float64}`

Default: `Vector{Float64}()`. Allows the user to pass a vector with the initial condition values desired in the simulation. If `initialize_simulation = true`, these values are used as a first guess and overwritten.

### `frequency_reference`

Default: `ReferenceBus`. Determines which frequency model is used for the network. Currently there are two options available:
- `ConstantFrequency` assumes that the network frequency is 1.0 per unit at all times.
- `ReferenceBus` will use the frequency state of a Dynamic Generator (rotor speed) or Dynamic Inverter (virtual speed) connected to the Reference Bus (defined in the Power Flow data) as the network frequency. If multiple devices are connected to such bus, the device with larger base power will be used as a reference. If a Voltage Source is connected to the Reference Bus, then a `ConstantFrequency` model will be used.

Future implementations will include a `CenterOfInertia` frequency model, and a `FrequencyDivider` model.

### `all_lines_dynamic::Bool`

Default: `false`. If `true` transforms all lines models (admittance matrix) into dynamic lines. See the [Network](component_models/network.md) page for more details.

### `all_branches_dynamic::Bool`

Default: `false`. If `true` transforms all branches (lines + transformers) into dynamic branches. See the [Network](component_models/network.md) page for more details.

### `system_to_file::Bool`

Default: `false`. Serializes the initialized system in the Folder Directory.

### `console_level::Logging`

Default: `Logging.Warn`. Sets the level of logging output to the console. Can be set to `Logging.Error`, `Logging.Warn`, `Logging.Info` or `Logging.Debug`.

### `file_level::Logging`

Default: `Logging.Debug`. Sets the level of logging output to file. Can be set to `Logging.Error`, `Logging.Warn`, `Logging.Info` or `Logging.Debug`

### `disable_timer_output::Bool`

Default: `false`. Allows the user to display timer information about the construction and initilization of the Simulation. If `true` disables the timer information.