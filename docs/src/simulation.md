# Simulation Structure

After constructing the System data from `PowerSystems.jl` with its dynamic components, a `Simulation` structure must be constructed. To do so, the arguments given are the following:

```julia
function Simulation
    ::SimulationModel
    system::PowerSystems.System
    simulation_folder::String
    tspan::NTuple{2, Float64}
    perturbations::Vector{<:Perturbation} = Vector{Perturbation}();
    kwargs...,
end
```

The Simulation will construct all the indexing for the states and will attempt to initialize the dynamical system in steady state based on the power flow data provided.

**Arguments:**
- `::SimulationModel` : Type of Simulation Model. `ResidualModel` or `MassMatrixModel`. See [Models Section](models.md) for more details
- `system::PowerSystems.System` : System data
- `simulation_folder::String` : Folder directory
- `tspan::NTuple{2, Float64}` : Time span for simulation
- `perturbations::Vector{<:Perturbation}` : Vector of Perturbations for the Simulation. Default: No Perturbations

**Accepted Keywords:**
- `initialize_simulation::Bool` : Runs the initialization routine. If false, simulation runs based on the operation point stored in System`
- `initial_conditions::Vector{Float64}` : Allows the user to pass a vector with the initial condition values desired in the simulation. If initialize_simulation = true, these values are used as a first guess and overwritten.
- `frequency_reference` : Default `ReferenceBus`. Determines which frequency model is used for the network. Currently there are two options available: 
    - `ConstantFrequency` assumes that the network frequency is 1.0 per unit at all times.
    - `ReferenceBus` will use the frequency state of a Dynamic Generator (rotor speed) or Dynamic Inverter (virtual speed) connected to the Reference Bus (defined in the Power Flow data) as the network frequency. If multiple devices are connected to such bus, the device with larger base power will be used as a reference. If a Voltage Source is connected to the Reference Bus, then a `ConstantFrequency` model will be used. 
    - Future implementations will include a `CenterOfInertia` frequency model, and a `FrequencyDivider` model
- `system_to_file::Bool` : Default `false`. Serializes the initialized system
- `console_level::Logging` : Default `Logging.Warn`. Sets the level of logging output to the console. Can be set to Logging.Error, Logging.Warn, Logging.Info or Logging.Debug
- `file_level::Logging` : Default `Logging.Debug`. Sets the level of logging output to file. Can be set to Logging.Error, Logging.Warn, Logging.Info or Logging.Debug
- `disable_timer_output::Bool` : Default `false`. Allows the user to display timer information about the construction and initilization of the Simulation.