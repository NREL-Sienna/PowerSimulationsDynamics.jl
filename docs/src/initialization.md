# Initialization Routine

Dynamic Simulations require a reasonable initial condition for the model. In most analysis,
power systems models are initialized at a stable equilibrium, which implies that:

```math
\begin{align}
0 = F(x, u, \eta)
\end{align}
```

Finding the solution of a large non-linear system is challenging and requires a reasonable
initial guess. In classical power systems literature, the routine to find equilibrium points
for the dynamic injection devices' components is well known and used in free and commercial
software (see [Power System Modelling and Scripting](https://www.springer.com/gp/book/9783642136689) page 224).
However, in the case of converter interface dynamic injection models, such routines are not documented.

This page shows the initialization details in `PowerSimulationsDynamics.jl`

## System-wide routine

The initialization routine starts from the solution of the power flow equations. For each
dynamic injection device `PowerSimulationsDynamics.jl` finds the solution of the systems of
non-linear equations for each dynamic component following the sequences described in the forthcoming
sections.

Once each device is individually initialized, the system-wide initial guess is used to solve the
system (1). In a first attempt at finding the solution, the tolerance is set to a stringent
tolerance. If the non-linear solver is unable to get a solution, it might usually reflect
small signal stability problems in the system. In a second attempt, the tolerances are relaxed.
If the solver succeeds, the simulation continues, but the user is warned.

```@raw html
<img src="./assets/sys_init.png" width="65%"/>
```

## Initialization of the Synchronous Machines

The initialization of Synchronous Machines is standard in power systems and follows the scheme
shown in the figure. Other internal variables are calculated recursively from the power flow
solution for the node on which the dynamic device isconnected. (Adapted from
[Power System Modelling and Scripting](https://www.springer.com/gp/book/9783642136689) Figure 9.2)

![init_machine](assets/synch_init.png)

## Initialization of the Inverters

Initializing the inverters follows the sequence shown in the figure and has been developed
to be compatible with the implementation of custom dynamic components. Given that the process
is less studied and standard than the one of the synchronous machine, this page contains more
detailed documentation of the process.

![init_machine](assets/inverter_init.png)

1. The first component to be initialized is the filter. Given that the filter is an RLC
   circuit connected to the grid, its currents and voltages need to match the results of the
   power flow. The initialization of the filter provides the values for the ``P`` and ``Q``
   used in the outer control and the ``V`` and ``I`` needed in the inner controls.
2. Based on the bus voltage in the system's reference frame ``V_r`` and the bus angle ``\theta``
   the PLL's can be initialized to obtain the angle and frequency estimates needed by the
   outer control.
3. The Outer Control calculates the internal angle ``\delta_{olc}`` required by the inner control
   to estimate the voltage and current phase difference.
4. The DC Source uses the power set-point consistent with the power outputs of the filter
   to initialize the ``V_{dc}`` set-points. This value is used in the inner control.
5. The inner control takes the phase angle ``\delta_{olc}`` and the ``V_{dc}`` to estimate the
   modulation values of the PWM converter.

**Note:** The initialization of an inverter through the proposed meta-model is actively under
development and subject to change. This page will maintain the latest version of the sequence.
