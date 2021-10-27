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

### Example

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

Another perturbation that can be modeled is a three phase fault at Bus 1 with impedance `r_f = 0.0001, x_f = 0.0` per unit, then the admittance of this new system is:

```julia
new_yb2 = [10000.0 - 20.0im  0.0 + 20.0im
           0.0 + 20.0im  0.0 - 20.0im]
```
Then, this perturbation ocurring at ``t = 1.0`` seconds can be included as:

```julia
ns2 = NetworkSwitch(1.0, new_yb2)
```

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