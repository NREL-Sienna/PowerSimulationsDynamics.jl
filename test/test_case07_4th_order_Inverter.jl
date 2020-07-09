"""
Case 7:
This case study a three bus system with 1 machine (One d- One q-: 4th order model), a VSM of 19 states and an infinite source.
The perturbation increase the reference power (analogy for mechanical power) of the machine from 0.6 to 0.8.
"""

##################################################
############### LOAD DATA ########################
##################################################

include(joinpath(dirname(@__FILE__), "data_tests/test07.jl"))

##################################################
############### SOLVE PROBLEM ####################
##################################################

#time span
tspan = (0.0, 20.0);
case_inv = collect(PSY.get_components(PSY.DynamicInverter, threebus_sys))[1]

#Define Fault using Callbacks
Pref_change = LITS.ControlReferenceChange(1.0, case_inv, LITS.P_ref_index, 1.2)

#Define Simulation Problem
sim = LITS.Simulation(threebus_sys, tspan, Pref_change)

small_sig = small_signal_analysis(sim)

#Solve problem in equilibrium
run_simulation!(sim, IDA(), dtmax = 0.02);

#Obtain data for angles
series = get_state_series(sim, ("generator-3-1", :θ_oc));

@test LinearAlgebra.norm(sim.x0_init - test07_x0_init) < 1e-6
@test sim.solution.retcode == :Success
@test small_sig.stable
