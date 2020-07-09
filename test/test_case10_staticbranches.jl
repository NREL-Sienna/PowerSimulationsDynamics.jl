"""
Case 8:
This case study a three bus system with 1 machine (One d- One q-: 4th order model), a VSM of 19 states and an infinite source. All lines are modeled as a static lines.
The perturbation trips two of the three circuits of line between buses 1 and 2, triplicating its impedance.
"""

##################################################
############### LOAD DATA ########################
##################################################

include(joinpath(dirname(@__FILE__), "data_tests/test08.jl"))

##################################################
############### SOLVE PROBLEM ####################
##################################################

#time span
tspan = (0.0, 40.0)

#Define Fault: Change of YBus
Ybus_change = LITS.ThreePhaseFault(
    1.0, #change at t = 1.0
    Ybus_fault,
) #New YBus

#Define Simulation Problem
sim = Simulation(
    threebus_sys, #system
    tspan, #time span
    Ybus_change, #Type of Fault
)

#Obtain small signal results for initial conditions
small_sig = small_signal_analysis(sim)

#Solve problem in equilibrium
run_simulation!(sim, IDA());

#Obtain data for voltages
series = get_voltagemag_series(sim, 2)

zoom = [
    (series[1][ix], series[2][ix])
    for (ix, s) in enumerate(series[1]) if (s > 0.90 && s < 1.6)
]

@test LinearAlgebra.norm(sim.x0_init - test10_x0_init) < 1e-6
@test sim.solution.retcode == :Success
@test small_sig.stable
