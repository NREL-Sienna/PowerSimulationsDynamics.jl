"""
Case 7:
This case study a three bus system with 1 machine (One d- One q-: 4th order model), a VSM of 19 states and an infinite source.
The perturbation increase the reference power (analogy for mechanical power) of the machine from 0.6 to 0.8.
"""

##################################################
############### LOAD DATA ########################
##################################################

##################################################
############### SOLVE PROBLEM ####################
##################################################

#time span
tspan = (0.0, 20.0);

x0_guess = [
    1.00, #V1_R
    1.00, #V2_R
    1.00, #V3_R
    0.0, #V1_I
    -0.01, #V2_I
    -0.01, #V3_I
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
    -0.1, #Iq_filter
    1.0, #eq_p
    0.47, #ed_p
    0.6, #δ
    1.0, #ω
    2.1, #Vf
    0.28, #Vr1
    -0.39, #Vr2,
    1.0,
] #Vm

#Define Fault using Callbacks
Pref_change = LITS.ControlReferenceChange(1.0, case7_gen, LITS.P_ref_index, 0.8)

#Define Simulation Problem
sim = LITS.Simulation(sys, tspan, Pref_change, initial_guess = x0_guess)

#Solve problem in equilibrium
run_simulation!(sim, IDA());

#Obtain data for angles
series = get_state_series(sim, ("Case7Gen", :δ));

@test sim.solution.retcode == :Success
