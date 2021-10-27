# Perturbations

Perturbations are used to alter the system from its steady state operation. If a Simulation is properly initialized, all states will remain fixed in their initial condition if no perturbation is applied to the system. The list of perturbations is as follows:

## NetworkSwitch

A `NetworkSwitch` allows to modify directly the admittance matrix, `Ybus`, used in the Simulation. This allows the user to perform branch modifications, three phase faults (with impedance larger than zero) or branch trips, as long as the new `Ybus` provided captures that perturbation.

The constructor is given by:

```raw
ns = NetworkSwitch(
    Time of Ocurrence,
    New Ybus
)
```

The time of ocurrence, defined as `Float64`, defines when the Network Switch will happen. This time should be inside the time span considered in the Simulation. The New Ybus, defined as a `Matrix{Float64}`, provides the new Ybus.

### Example 1: Circuit Disconnection

Consider a two bus system connected via a double circuit line, on which each circuit has parameters, `r = 0.0, x = 0.1, b = 0.0` per unit, then the admittance matrix of the original system is given by:

```julia
yb = [0.0 - 20.0im 0.0 + 20.0im
      0.0 + 20.0im 0.0 - 20.0im]
```

Triping one circuit can be modeled by doubling the impedance, i.e., dividing by 2 the admittance:

```julia
new_yb = [0.0 - 10.0im 0.0 + 10.0im
          0.0 + 10.0im 0.0 - 10.0im]
```

Then, this perturbation ocurring at ``t = 1.0`` seconds can be included as:
```julia
ns1 = NetworkSwitch(1.0, new_yb)
```

### Example 2: Three Phase Fault

Another perturbation that can be modeled is a three phase fault at Bus 1 with impedance `r_f = 0.0001, x_f = 0.0` per unit, then the admittance of this new system is:

```julia
new_yb2 = [10000.0 - 20.0im  0.0 + 20.0im
           0.0 + 20.0im  0.0 - 20.0im]
```
Then, this perturbation ocurring at ``t = 1.0`` seconds can be included as:

```julia
ns2 = NetworkSwitch(1.0, new_yb2)
```

Now, consider that the fault is cleared at ``t = 1.05`` seconds by disconnecting the Circuit 2 of the line. This can be modeled with the single circuit admittance matrix:

```julia
new_yb3 = [0.0 - 10.0im 0.0 + 10.0im
           0.0 + 10.0im 0.0 - 10.0im]
```

and the perturbation as:

```julia
ns3 = NetworkSwitch(1.05, new_yb3)
```

Then, the entire perturbation for the Simulation can be included in a vector of perturbations as:

```julia
three_fault = [ns2, ns3]
```

that can be passed as a perturbation argument in the Simulation construction.

## BranchTrip

A `BranchTrip` completely disconnects a branch from the system. Currently there is only support for static branches disconnection. Future releases will provide support for a Dynamic Line disconnection. The constructor is as follows:

```raw
ns = BranchTrip(
    Time of Ocurrence,
    Type of Line,
    Name of Line,
)
```

The time of ocurrence, defined as `Float64`, defines when the Branch Trip will happen. The supported type of branches can be `Line` and `Transformer2W`. `DynamicLine` is currently not supported. The Name of Line, defined as `String`, determines which branch will be disconnected

### Example

Consider the following 2 bus system defined by:

```julia
buses = [
    Bus(1, "nodeA", "REF", 0, 1.0, (min = 0.9, max = 1.05), 230, nothing, nothing),
    Bus(2, "nodeB", "PV", 0, 1.0, (min = 0.9, max = 1.05), 230, nothing, nothing),
]

line1 = Line(
        "Circuit1",
        true,
        0.0,
        0.0,
        Arc(from = buses[1], to = buses[2]),
        0.00,
        0.1,
        (from = 0.0, to = 0.0),
        2.0,
        (min = -0.7, max = 0.7),
    )
line2 = Line(
        "Circuit2",
        true,
        0.0,
        0.0,
        Arc(from = buses[1], to = buses[2]),
        0.0,
        0.1,
        (from = 0.0, to = 0.0),
        2.0,
        (min = -0.7, max = 0.7),
    )

sys = System(100.0, buses, [], [], [line1, line2])
```

A Branch Trip of Circuit 2 at time ``t = 1.0`` seconds, can be implemented as:

```julia
b_trip = BranchTrip(1.0, Line, "Circuit2")
```

**Note:** Islanding is currently not supported in `PowerSimulationsDynamics.jl`. If a `BranchTrip` isolates a generation unit, the system may diverge due to the isolated generator.

## GeneratorTrip

A `GeneratorTrip` allows to disconnect a Dynamic Generation unit from the system at a specified time. The constructor is the following

```raw
g_trip = GeneratorTrip(
    Time of Ocurrence,
    Dynamic Device,
)
```

The time of ocurrence, defined as `Float64`, defines when the Generator Trip will happen. This time should be inside the time span considered in the Simulation. The Dynamic Device, defined as `DynamicInjection`, must be taken from the system to identify which device will be disconnected from the system.

### Example

onsider that you have a generator at bus 102, named `"generator-102-1"` in your system called `sys`. The constructor to trip it from the system is:

```julia
g = get_component(DynamicGenerator, sys, "generator-102-1")
g_trip = GeneratorTrip(1.0, g)
```

## ControlReferenceChange

A `ControlReferenceChange` allows to change the reference setpoint provided by a generator/inverter. The constructor is the following:

```raw
crc = ControlReferenceChange(
    Time of Ocurrence,
    Dynamic Device,
    Signal,
    New Reference Value
)
```

The time of ocurrence, defined as `Float64`, defines when the Control Reference Change will happen. This time should be inside the time span considered in the Simulation. The Dynamic Device, defined as `DynamicInjection`, must be taken from the system to identify which device will modify its setpoint. The Signal, defined as `Symbol`, determines which reference setpoint will be modified. The accepted signals are:
- `:P_ref`: Modifies the active power reference setpoint.
- `:V_ref`: Modifies the voltage magnitude reference setpoint.
- `:Q_ref`: Modifies the reactive power reference setpoint (if used).
- `:ω_ref`: Modifies the frequency setpoint.

The New Reference Value, defined as `Float64`, updates the setpoint.

### Example

Consider that you have a generator at bus 102, named `"generator-102-1"` in your system called `sys`. The constructor to change is active power reference to `0.5` is:

```julia
g = get_component(DynamicGenerator, sys, "generator-102-1")
crc = ControlReferenceChange(1.0, g, :P_ref, 0.5)
```

## LoadChange

A `LoadChange` allows to change the active or reactive power setpoint from a load. The constructor is the following:

```raw
l_change = LoadChange(
    Time of Ocurrence,
    Load Device,
    Signal,
    New Reference Value
)
```

The time of ocurrence, defined as `Float64`, defines when the Load Change will happen. This time should be inside the time span considered in the Simulation. The Load Device, defined as `ElectricLoad`, must be taken from the system to identify which load will modify its setpoint. The Signal, defined as `Symbol`, determines which reference setpoint will be modified. The accepted signals are:
- `:P_ref`: Modifies the active power reference setpoint.
- `:Q_ref`: Modifies the reactive power reference setpoint. 

The New Reference Value, defined as `Float64`, updates the setpoint.

### Example

Consider that you have a load at bus 103, named `"load-103-1"` in your system called `sys`. The constructor to change is active power reference to `0.8` per unit at ``t = 1.0`` seconds is:

```julia
l_device = get_component(ElectricLoad, sys, "load-103-1")
l_change = LoadChange(1.0, l_device, :P_ref, 0.8)
```

## LoadTrip

A `LoadTrip` allows the user to disconnect a load from the system. The constructor is the following:

```raw
l_trip = LoadTrip(
    Time of Ocurrence,
    Load Device,
)
```

The time of ocurrence, defined as `Float64`, defines when the Load Change will happen. This time should be inside the time span considered in the Simulation. The Load Device, defined as `ElectricLoad`, must be taken from the system to identify which load will be disconnected from the system.

### Example

Consider that you have a load at bus 103, named `"load-103-1"` in your system called `sys`. The constructor to disconnect such load at ``t = 1.0`` seconds is:

```julia
l_device = get_component(ElectricLoad, sys, "load-103-1")
l_trip = LoadTrip(1.0, l_device)
```

## SourceBusVoltageChange

A `SourceBusVoltageChange` allows to change the reference setpoint provided by a voltage source. The constructor is the following:

```raw
sbvc = SourceBusVoltageChange(
    Time of Ocurrence,
    Source Device,
    Signal Index,
    New Reference Value
)
```

The time of ocurrence, defined as `Float64`, defines when the Source Voltage Change will happen. This time should be inside the time span considered in the Simulation. The Source Device, defined as `Source`, must be taken from the system to identify which load will modify its setpoint. The Signal Index, defined as `Int`, determines which reference setpoint will be modified. The accepted signals are:
- `1`: Modifies the internal voltage magnitude reference setpoint.
- `2`: Modifies the internal voltage angle reference setpoint.

The New Reference Value, defined as `Float64`, updates the setpoint.

### Example

Consider that you have a voltage source at bus 101, named `"source-101-1"` in your system called `sys`. The constructor to change is voltage magnitude reference to `1.02` per unit at ``t = 1.0`` seconds is:

```julia
s_device = get_component(Source, sys, "source-101-1")
s_change = SourceBusVoltageChange(1.0, s_device, 1, 1.02)
```