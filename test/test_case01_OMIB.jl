"""
Case 1:
This case study defines a classical machine against an infinite bus. The fault
drop a circuit on the (double circuit) line connecting the two buses, doubling its impedance
"""

##################################################
############### LOAD DATA ########################
##################################################

include(joinpath(dirname(@__FILE__), "data_tests/test01.jl"))

##################################################
############### SOLVE PROBLEM ####################
##################################################
#Define Fault: Change of YBus
Ybus_change = ThreePhaseFault(
    1.0, #change at t = 1.0
    Ybus_fault,
) #New YBus

#Define Simulation Problem
sim = Simulation(
    omib_sys, #system
    (0.0, 30.0), #time span
    Ybus_change,
) #Type of Fault

#Obtain small signal results for initial conditions
small_sig = small_signal_analysis(sim)

#Solve problem in equilibrium
run_simulation!(sim, IDA(), dtmax = 0.02);

#Obtain data for angles
series = get_state_series(sim, ("generator-2-1", :δ));
series2 = get_voltagemag_series(sim, 2)
LITS.print_init_states(sim)

@test LinearAlgebra.norm(sim.x0_init - test01_x0_init) < 1e-6
@test sim.solution.retcode == :Success
@test small_sig.stable
