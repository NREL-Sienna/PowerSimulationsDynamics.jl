"""
Case 6:
This case study a 19-state virtual synchronous machine against an infinite bus located at bus 1, with VSM located at bus 2.
The perturbation increase the reference power (analogy for mechanical power) from 0.5 to 0.7.
"""

##################################################
############### LOAD DATA ########################
##################################################

############### Data Network ########################

nodes_DAIB = nodes_DArco_IB()

branch_DAIB = branches_DArco_IB(nodes_DAIB)

############### Data devices ########################

inf_gen_DAIB = inf_gen_1_pu(nodes_DAIB)

############### Inverter Data ########################

Darco_Inverter = inv_DAIB(nodes_DAIB)

######################### Dynamical System ########################

#Create system with BasePower = 100 MVA and nominal frequency 50 Hz.
sys = system_DAIB(nodes_DAIB, branch_DAIB, [inf_gen_DAIB], [Darco_Inverter])

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
    0.0, #δω_vsm
    0.2, #δθ_vsm
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
    0.1, #δθ_pll
    0.5, #id_cv
    0.0, #iq_cv
    0.95, #vod
    -0.1, #voq
    0.49, #iod
    -0.1,
] #ioq

#Define Fault using Callbacks
Pref_change = LITS.ControlReferenceChange(1.0, Darco_Inverter, LITS.P_ref_index, 0.7)

#Define Simulation Problem
sim = LITS.Simulation(sys, tspan, Pref_change, initial_guess = x0_guess)

#Solve problem in equilibrium
run_simulation!(sim, Sundials.IDA());

#Obtain data for angles
series = get_state_series(sim, ("DARCO", :δω_vsm))

@test sim.solution.retcode == :Success
