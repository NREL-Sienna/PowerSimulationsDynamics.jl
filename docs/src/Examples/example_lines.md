# Tutorial: Dynamic Lines

This tutorial will introduce an example of considering dynamic lines in `LITS`.
Note that this tutorial is for `LITS 0.4.0`. 

This tutorial presents a simulation of a three-bus system, with an infinite bus (represented as a voltage source behind an impedance) at bus 1, a one d- one q- machine on bus 2 and an inverter of 19 states, as a virtual synchronous machine at bus 3. The perturbation will be the trip of two of the three circuits (triplicating its resistance and impedance) of the line that connects bus 1 and bus 3. This case also consider a dynamic line model for connection between buses 2 and 3.

It is recommended to check `Tutorial 1: OMIB` first, since that includes more details and explanations on all definitions and functions.

This tutorial can be found on [LITS/Examples](https://github.com/Energy-MAC/LITS-Examples) repository.

## Step 1: Package Initialization

```julia
using LITS
using PowerSystems
using Sundials
const PSY = PowerSystems
```

## Step 2: Data creation

We will load the system already defined in Tutorial 0:

```julia
threebus_sys = System("threebus_sys.json")
```

In addition, we will create a new copy of the system on which we will simulate the same case, but will consider dynamic lines:
```julia
threebus_sys_dyn = deepcopy(threebus_sys)
```

## Step 3: Build the simulation and initializing the problem

First, we construct the perturbation, by properly computing the new Ybus:

```julia
#Create Ybus_Fault
fault_branches = deepcopy(collect(get_components(Branch, threebus_sys)))
for br in fault_branches
    if get_name(br) == "1"
        br.r = 3 * br.r
        br.x = 3 * br.x
        b_new = (from = br.b.from / 3, to = br.b.to / 3)
        br.b = b_new
    end
end
Ybus_fault = Ybus(fault_branches, get_components(Bus, threebus_sys))[:, :];
#Define Fault: Change of YBus
Ybus_change = LITS.ThreePhaseFault(
    1.0, #change at t = 1.0
    Ybus_fault, #New YBus
) 
```

Now, we construct the simulation:

```julia
#Time span of our simulation
tspan = (0.0, 30.0)

#Define Simulation
sim = Simulation(
    pwd(), #folder to output results
    threebus_sys, #system
    tspan, #time span
    Ybus_change, #Type of perturbation
)
```

We can obtain the initial conditions as:
```julia
#Will print the initial states. It also give the symbols used to describe those states.
print_device_states(sim)
#Will export a dictionary with the initial condition values to explore
x0_init = get_dict_init_states(sim)
```


## Step 4: Run the simulation

```julia
#Run the simulation
run_simulation!(sim, #simulation structure
                IDA(), #Sundials DAE Solver
                dtmax = 0.02, #Maximum step size
)
```

## Step 5: Explore the solution

```julia
using Plots
volt = get_voltagemag_series(sim, 102)
plot(volt, xlabel="time", ylabel="Voltage [pu]", label="V_2")

zoom = [(volt[1][ix], volt[2][ix]) for (ix, s) in enumerate(volt[1]) if (s > 0.90 && s < 1.6)]
plot(zoom, xlabel="time", ylabel="Voltage [pu]", label="V_2")
```

```@raw html
<img src="../../assets/voltage_lines.png" width="75%"/>
``` ⠀

```@raw html
<img src="../../assets/voltage_zoom_lines.png" width="75%"/>
``` ⠀
