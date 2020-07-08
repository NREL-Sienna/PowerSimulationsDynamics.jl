using LITS
using Sundials

"""
Case 6:
This case study a 19-state virtual synchronous machine against an infinite bus located at bus 1, with VSM located at bus 2.
The perturbation increase the reference power (analogy for mechanical power) from 0.5 to 0.7.
"""

##################################################
############### LOAD DATA ########################
##################################################

include(joinpath(dirname(@__FILE__), "data_tests/test06.jl"))

##################################################
############### SOLVE PROBLEM ####################
##################################################

#time span
tspan = (0.0, 4.0);

#Define Fault using Callbacks
Pref_change = LITS.ControlReferenceChange(1.0, case_inv, LITS.P_ref_index, 0.7)

#Define Simulation Problem
sim = LITS.Simulation(omib_sys, tspan, Pref_change)

small_sig = small_signal_analysis(sim)

#Solve problem in equilibrium
run_simulation!(sim, Sundials.IDA());

#Obtain data for angles
series = get_state_series(sim, ("generator-2-1", :ω_oc))

print_device_states(sim)

@test LinearAlgebra.norm(sim.x0_init - test06_x0_init) < 1e-6
@test sim.solution.retcode == :Success
@test small_sig.stable
