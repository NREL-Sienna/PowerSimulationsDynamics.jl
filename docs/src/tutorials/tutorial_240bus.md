#[PSSE 240 Bus Case system with Renewables](https://www.nrel.gov/grid/test-case-repository.html)

**Originally Contributed by**: José Daniel Lara

# Introduction

This tutorial will introduce the industry models of Renewable Energy the comparisons between DiffEq Integration techniques for comparison. We show the uses of Sundials and OrdinaryDiffEq to obtain the transient response of a system to a perturbation.

```@repl sys_240bus
using PowerSimulationsDynamics
using PowerSystemCaseBuilder
using PowerSystems
using Sundials
using Plots
using OrdinaryDiffEq
const PSY = PowerSystems
```

## Load the system

```@repl sys_240bus
file_dir = joinpath(pkgdir(PowerSimulationsDynamics), "test", "data_tests")
sys = System(joinpath(file_dir, "WECC_240_dynamic.json"); runchecks = false)
# Transform loads to constant impedance
include(joinpath(pkgdir(PowerSimulationsDynamics), "test/utils/get_results.jl"))
for l in get_components(PSY.StandardLoad, sys)
    transform_load_to_constant_impedance(l)
end
```

## Build the simulation and initialize the problem

The next step is to create the simulation structure. This will create the indexing of our system that will be used to formulate the differential-algebraic system of equations. To do so, it is required to specify the perturbation that will occur in the system. In this case, we will use a ResidualModel formulation, for more details about the formulation checkout the [Models Section](https://nrel-sienna.github.io/PowerSimulationsDynamics.jl/stable/models/) in `PowerSimulationsDynamics.jl` documentation.

```@repl sys_240bus
using Logging
sim_ida = Simulation(
    ResidualModel,
    sys, #system
    pwd(),
    (0.0, 20.0), #time span
    BranchTrip(1.0, Line, "CORONADO    -1101-PALOVRDE    -1401-i_10");
    console_level = Logging.Info,
)
```

## Run the simulation using Sundials

We will now run the simulation using Sundials.jl solver IDA() by specifying the maximum dt we want for the simulation. In our experience with this solver, solution times are faster when supplying information about the maximum time step than the tolerances as we can see in the example

```@repl sys_240bus
execute!(sim_ida, IDA(), dtmax = 0.01)
```

## Read the results and plot a system variable

After the simulation is completed, we can extract the results and make plots as desired. In this case, we will plot the voltage magnitude at the bus at which the line was connected.

```@repl sys_240bus
res_ida = read_results(sim_ida)
v1101_ida = get_voltage_magnitude_series(res_ida, 1101);
plot(v1101_ida);
```

![plot](figs/v1101_ida.svg)

## Run the simulation using Rodas4()

In this case, we will use a MassMatrixModel formulation, for more details about the formulation checkout the [Models Section](https://nrel-sienna.github.io/PowerSimulationsDynamics.jl/stable/models/) in `PowerSimulationsDynamics.jl` documentation

```@repl sys_240bus
sim_rodas = Simulation(
    MassMatrixModel,
    sys, #system
    pwd(),
    (0.0, 20.0), #time span
    BranchTrip(1.0, Line, "CORONADO    -1101-PALOVRDE    -1401-i_10");
    console_level = Logging.Info,
)
```

We will now run the simulation using OrdinaryDiffEq.jl solver Rodas4() by specifying the tolerance we want for the simulation. In our experience with this solver, solution times are faster when supplying information about the atol and rtol values as we can see in the example. The solver will also work with a specified dtmax but take a significantly longer time to solve. When using OrdinaryDiffEq.jl solvers always pass the option `initializealg = NoInit()` to avoid unnecessary re-initialization of the algebraic equations.

```@repl sys_240bus
execute!(
    sim_rodas,
    Rodas4(),
    saveat = 0.01,
    atol = 1e-10,
    rtol = 1e-10,
    initializealg = NoInit(),
)
```

## Read the results

After the simulation is completed, we can extract the results and make plots as desired. In this case, we will plot the voltage magnitude at the bus at which the line was connected.

```@repl sys_240bus
res_rodas = read_results(sim_rodas)
```

## Compare the results

After the simulation is completed, we can extract the results and make plots as desired. In this case, we will plot the voltage magnitude at the bus at which the line was connected. For both of the solution techniques.

```@repl sys_240bus
v1101 = get_voltage_magnitude_series(res_rodas, 1101);
plot(v1101, label = "RODAS4");
plot!(v1101_ida, label = "IDA");
```

![plot](figs/v1101_comparison.svg)