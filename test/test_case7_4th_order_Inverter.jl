"""
Case 7:
This case study a three bus system with 1 machine (One d- One q-: 4th order model), a VSM of 19 states and an infinite source.
The perturbation increase the reference power (analogy for mechanical power) of the machine from 0.6 to 0.8.
"""

##################################################
############### LOAD DATA ########################
##################################################

############### Data Network ########################

nodes_case7 = nodes_3bus()

branch_case7 = branches_3lines(nodes_case7)

#Trip of Line 1.
branch_case7_fault = branches_3lines_fault(nodes_case7)

loads_case7 = loads_3bus_case7(nodes_case7)

############### Data devices ########################

inf_gen_case7 = inf_gen_1_pu(nodes_case7)

### Case 7 Generators ###

case7_gen = dyn_gen_case7(nodes_case7)

############### Inverter Data ########################

case7_inv = inv_case78(nodes_case7)

######################### Dynamical System ########################

#Create system with BasePower = 100 MVA and nominal frequency 50 Hz.
sys = system_50Hz(
    nodes_case7,
    branch_case7,
    loads_case7,
    [inf_gen_case7],
    [case7_inv],
    [case7_gen],
)

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
    -0.1, #ioq
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
