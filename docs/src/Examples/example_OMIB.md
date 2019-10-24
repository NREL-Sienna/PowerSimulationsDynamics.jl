# Tutorial: One Machine against Infinite Bus (OMIB)

This tutorial will introduce you to the functionality of `LITS` for running Power System Simulations.
Note that this tutorial is for `LITS 0.2.0`. Future versions will have dedicated functions to find an equilibrium point and a proper functions for running perturbations without coding directly the callbacks.

This tutorial presents a simulation of a two-bus system, with an infinite bus (represented as a voltage source behind an impedance) at bus 1 and a classic machine on bus 2. The perturbation will be the trip of one of the two circuits (doubling its resistance and impedance) of the line that connects both buses.

This tutorial can be found on [LITS/Examples](https://github.com/Energy-MAC/LITS-Examples) repository.

## Step 1: Package Initialization

The first step consists in initialize all packages that will be used to run the simulation. All the necessary packages are listed:

```julia
using LITS
using PowerSystems
using NLsolve
using DiffEqBase
using Sundials
const PSY = PowerSystems
```

`LITS` and `PowerSystems` are used to properly define the data structure, while `NLsolve`, `DiffEqBase` and `Sundials` are used to formulate and solve the problem. Finally we call use can call `PowerSystems` functions using the `PSY` abbreviation.

## Step 2: Data creation

Next we need to define the different elements required to run a simulation. To run a simulation, it is required to define a `DynamicSystem` that requires the following components:

- Vector of `PSY.Bus` elements, that define all the buses in the network.
- Vector of `PSY.Branch` elements, that define all the branches elements (that connect two buses) in the network.
- Vector of `DynInjection` elements, that define all the devices connected to buses that can inject (or withdraw) current, while also defining differential equations to model its dynamics.
- Vector of `PSY.Injection` elements, that define all the devices connected to buses that can inject (or withdraw) current, without defining any differential equation.
- The base of power used to define per unit values, in MVA as a `Float64` value.
- The base frequency used in the system, in Hz as a `Float64` value.
- (Optional) Vector of `DynBranch` elements, that can be used to model lines with differential equations.

To start we will define the data structures for the network.

### Buses and Branches

As mentioned earlier, we require to create a `Vector` of `PSY.Bus` to define the buses in the network. Currently, some of the parameters are not used in `LITS`, but will be used once the initialization procedure is implemented.

```julia
#Define the vector of buses
nodes_case1 = [PSY.Bus(1 , #number
                   "Bus 1", #Name
                   "REF" , #BusType (REF, PV, PQ)
                   0, #Angle in radians
                   1.05, #Voltage in pu
                   (min=0.94, max=1.06), #Voltage limits in pu
                   69), #Base voltage in kV
                   PSY.Bus(2 , "Bus 2"  , "PV" , 0 , 1.0 , (min=0.94, max=1.06), 69)];
```

Note that two buses are defined in the vector `nodes_case1`. Similarly, to define the branches (that also has some parameters that are currently not used):

```julia
#Define the vector of branches
branch_case1 = [PSY.Line("Line1", #name
                     true, #available
                     0.0, #active power flow initial condition (from-to)
                     0.0, #reactive power flow initial condition (from-to)
                     Arc(from=nodes_case1[1], to=nodes_case1[2]), #Connection between buses
                     0.01, #resistance in pu
                     0.05, #reactance in pu
                     (from=0.0, to=0.0), #susceptance in pu
                     18.046, #rate in MW
                     1.04)];  #angle limits (-min and max)
```

Since we are interested in creating a fault that trips one of the two circuits of the line, we will create an additional `Vector` of branches with doubled impedance:

```julia
#Define the vector of branches under the fault
branch_case1_fault = [PSY.Line("Line1", #name
                           true, #available
                           0.0, #active power flow initial condition (from-to)
                           0.0, #reactive power flow initial condition (from-to)
                           Arc(from=nodes_case1[1], to=nodes_case1[2]), #Connection between buses
                           0.02, #resistance in pu
                           0.1, #reactance in pu
                           (from=0.0, to=0.0), #susceptance in pu
                           18.046, #rate in MW
                           1.04)];  #angle limits (-min and max)
```

Note that the resistance and reactance is doubled in comparison to the system without fault.

### Injection devices

Secondly, we will define devices that can inject/withdraw electric current directly without defining differential equations. In this case we include a load and the voltage source that model the infinite bus.

```julia
loads_case1 = [PowerLoad("Bus1", #name
                         true, #available
                         nodes_case1[2], #bus location
                         PowerSystems.ConstantPower, #type
                         0.3, #Active Power pu
                         0.01, #Reactive power pu
                         0.3, #Max Active Power pu
                         0.01)]; #Max Reactive Power pu

inf_gen_case1 = StaticSource(1, #number
                :InfBus, #name
                nodes_case1[1], #bus
                1.05, #VR real part of voltage source
                0.0, #VI imaginary part of voltage source
                0.000001); #Xth
```

Note that loads are assumed as constant power for power flow purposes, but for dynamic simulations are converted to impedance loads assuming nominal voltage equals to 1 pu.

### Dynamic Injection devices

Third, we define the `Vector` of `DynInjection` elements. In this case, we require to define a generator located in bus 2. For that purpose, we need to define its machine, shaft, automatic voltage regulator (AVR), turbine governor (TG) and power system stabilizer (PSS):

```julia
### Machine ###
case1_machine = BaseMachine(0.0, #R
                            0.2995, #Xd_p
                            0.7087, #eq_p
                            100.0)  #MVABase

######## Shaft Data #########

### Shaft for Case 1 ###
case1_shaft = SingleMass(3.148, #H
                     2.0) #D

########  AVR Data #########
case1_avr = AVRFixed(0.0) #Vf not applicable in Classic Machines

######## TG Data #########
### No TG ###
case1234_no_tg = TGFixed(1.0) #eff

######## PSS Data #########
### No PSS ###
cases_no_pss = PSSFixed(0.0) #No PSS without AVR

### Constructing the Generator ###
case1_gen = DynGenerator(1, #Number
                      :Case1Gen,
                      nodes_case1[2], #bus
                      1.0, # ω_ref,
                      1.0, #V_ref (only used in AVR)
                      0.5, #P_ref
                      case1_machine, #machine
                      case1_shaft, #shaft
                      case1_avr, #avr
                      case1234_no_tg, #tg
                      cases_no_pss) #pss
```

Note that a generator is defined by its 5 components, while also defining its reference for frequency, voltage and power.

### Defining the Dynamic System

Finally, with all the components properly constructed we define the dynamic system:

```julia
case1_DynSystem = DynamicSystem(nodes_case1, #Vector of Buses
                              branch_case1, #Vector of Branches
                              [case1_gen], #Vector of Dynamic Injections
                              vcat(inf_gen_case1,loads_case1), #Vector of Injections
                              100.0, #MVA Base
                              60.0) #Freq. Base
```

## Step 3: Initializing the problem

The next step consists in finding an initial condition for the states. But first, we will explore some of the characteristics of our Dynamic System. All information (a ton) can be observed using `dump(case1_DynSystem)`. The following methods can be used to return some information:

- `case1_DynSystem.buses`: Return the vector of buses of the dynamic system.
- `case1_DynSystem.branches`: Return the vector of branches of the dynamic system.
- `case1_DynSystem.dyn_injections`: Return the vector of dynamic injections of the dynamic system.
- `case1_DynSystem.injections`: Return the vector of dynamic injections of the dynamic system.
- `case1_DynSystem.DAE_vector`: Return the vector of booleans of the dynamic system. Returns `false` or 0 for states that are algebraic and `true` or 1 for states that have derivative defined (differential states). The arrangement will put first the real part of the voltage buses and next the imaginary part. After that the differential states are defined.
- `case1_DynSystem.global_state_index`: Return an array of dictionaries that have the order of the states in the entire vector state.

To initialize the problem we need to define an initial guess of the states:

```julia
#Initialize variables
dx0 = zeros(LITS.get_total_rows(case1_DynSystem)) #Define a vector of zeros for the derivative
x0 = [1.05, #VR_1
      1.0, #VR_2
      0.0, #VI_1
      0.01, #VI_2
      0.2, #δ
      1.0] #ω
tspan = (0.0, 30.0);
```

We will use `NLsolve` to find the initial condition of the system:

```julia
inif! = (out,x) -> LITS.system_model!(out, #output of the function
                                      dx0, #derivatives equal to zero
                                      x, #states
                                       ([0.0],case1_DynSystem), #Parameters: [0.0] is not used
                                        0.0) #time equals to zero.
sys_solve = nlsolve(inif!, x0) #Solve using initial guess x0
x0_init = sys_solve.zero
```

## Step 4: Build the Simulation

Next we will construct the simulation that we are interested to run. But first, we define the pertubation we are interested in model. `LITS` have two perturbations already implemented, that are a change in the mechanical power `P_ref` and a change on the admittance matrix `Y_bus` of the system. In this case we define a change in the admittance matrix:

```julia
#Compute Y_bus after fault
Ybus_fault = PSY.Ybus(branch_case1_fault, nodes_case1)[:,:] #Obtain Ybus for fault system

#Define Fault using Callbacks
cb = DiffEqBase.DiscreteCallback(LITS.change_t_one, #Change occurs at t=1
                                 LITS.Y_change!) #Callback will change the Y_bus.
```

Now we define the simulation structure:

```julia
#Define Simulation Problem
sim = DynamicSimulation(case1_DynSystem, #Dynamic System
                        tspan, #Time span to simulate
                        Ybus_fault, #Parameter that will be changed in the fault
                        cb, #Callback
                        x0_init) #Initial condition
```

Finally, to run the simulation:
```julia
#Solve problem in equilibrium
run_simulation!(sim, #simulation structure
                IDA(), #Sundials DAE Solver
                dtmax=0.02); #Arguments: Maximum timestep allowed
```

## Step 5: Exploring the solution

After running the simulation, our simulation structure `sim` will have the solution. For that `sim.solution` can be used to explore the solution structure. In this case `sim.solution.t` returns the vector of time, while `sim.solution.u` return the array of states. In addition, `LITS` have two functions to obtain different states of the solution:

- `get_state_series(sim, (:Case1Gen, :δ))`: can be used to obtain the solution as a tuple of time and the required state. In this case, we are obtaining the rotor angle `:δ` of the generator named `:Case1Gen`.
- `get_voltagemag_series(sim, 2)`: can be used to obtain the voltage magnitude as a tuple of time and voltage. In this case, we are obtaining the voltage magnitude at bus 2 (where the generator is located).

```julia
using Plots
angle = get_state_series(sim, (:Case1Gen, :δ))
plot(angle, xlabel="time", ylabel="rotor angle [rad]", label="rotor angle")

volt = get_voltagemag_series(sim, 2)
plot(volt, xlabel="time", ylabel="Voltage [pu]", label="V_2")
```

```@raw html
<img src="../../assets/rotor_angle_OMIB.png" width="75%"/>
``` ⠀

```@raw html
<img src="../../assets/voltage_OMIB.png" width="75%"/>
``` ⠀
