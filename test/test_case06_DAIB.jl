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

include(joinpath("test", "data_tests/test06.jl"))

##################################################
############### SOLVE PROBLEM ####################
##################################################

#time span
tspan = (0.0, 4.0);

#Initial guess
x0_guess = [
    1.00, #V1_R
    1.0648, #V2_R
    0.0, #V1_I
    0.001, #V2_I
    1.0, #ω_oc
    0.2, #θ_oc
    0.025, #qm
    0.0015, #ξ_d
    -0.07, #ξ_q
    0.05, #γ_d
    -0.001, #γ_q
    0.95, #ϕ_d
    -0.10, #ϕ_q
    1.004, #vpll_d
    0.0, #vpll_q
    0.0, #ε_pll
    0.1, #θ_pll
    0.5, #id_cv
    0.0, #iq_cv
    0.95, #Vd_filter
    -0.1, #Vq_filter
    0.49, #Id_filter
    -0.1,
] #Iq_filter

#Define Fault using Callbacks
Pref_change = LITS.ControlReferenceChange(1.0, darco_inverter, LITS.P_ref_index, 0.7)

#Define Simulation Problem
sim = LITS.Simulation(omib_sys, tspan, Pref_change, initial_guess = x0_guess)

#Solve problem in equilibrium
run_simulation!(sim, Sundials.IDA());

#Obtain data for angles
series = get_state_series(sim, ("DARCO", :ω_oc))

@test sim.solution.retcode == :Success
