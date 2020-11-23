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

![init_system][assets/init.png]

## Initialization of the Synchronous Machines

The initialization of Synchronous Machines is standard in power systems and follows the scheme
shown in the figure. Other internal variables are calculated recursively from the power flow
solution for the node on which the dynamic device isconnected. (The figure is an adapatation
[Power System Modelling and Scripting](https://www.springer.com/gp/book/9783642136689) Figure 9.2)

![init_machine](assets/synch_init.png)

## Initialization of the Inverters



![init_machine](assets/inverter_init.png)
