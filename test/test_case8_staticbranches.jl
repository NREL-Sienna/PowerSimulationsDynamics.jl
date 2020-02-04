"""
Case 8:
This case study a three bus system with 1 machine (One d- One q-: 4th order model), a VSM of 19 states and an infinite source. All lines are modeled as a static lines.
The perturbation trips two of the three circuits of line between buses 1 and 2, triplicating its impedance.
"""

##################################################
############### LOAD DATA ########################
##################################################

############### Data Network ########################

nodes_case8 = nodes_3bus()

branch_case8 = branches_3lines_case8(nodes_case8)

#Two of three lines lost in Line 1
branch_case8_fault = branches_3lines_case8_fault(nodes_case8)

loads_case8 = loads_3bus_case8(nodes_case8)

############### Data devices ########################

inf_gen_case8 = inf_gen_1_pu(nodes_case8)

### Case 8 Generators ###

case8_gen = dyn_gen_case8(nodes_case8)

############### Inverter Data ########################

case8_inv = inv_case78(nodes_case8)

######################### Dynamical System ########################

#Create system with BasePower = 100 MVA and nominal frequency 50 Hz.
sys = system_50Hz(
    nodes_case8,
    branch_case8,
    loads_case8,
    [inf_gen_case8],
    [case8_inv],
    [case8_gen],
)

##################################################
############### SOLVE PROBLEM ####################
##################################################

#Compute Y_bus after fault
Ybus_fault = get_admittance_matrix(nodes_case8, branch_case8_fault)

#time span
tspan = (0.0, 20.0)

#initial_guess
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

#Define Fault: Change of YBus
Ybus_change = LITS.ThreePhaseFault(
    1.0, #change at t = 1.0
    Ybus_fault,
) #New YBus

#Define Simulation Problem
sim = Simulation(
    sys, #system
    tspan, #time span
    Ybus_change, #Type of Fault
    initial_guess = x0_guess,
) #initial guess

small_sig = small_signal_analysis(sim)

#Solve problem in equilibrium
run_simulation!(sim, IDA());

#Obtain data for voltages
series = get_voltagemag_series(sim, 2)

zoom = [
    (series[1][ix], series[2][ix])
    for (ix, s) in enumerate(series[1]) if (s > 0.90 && s < 1.6)
]

@test sim.solution.retcode == :Success
@test small_sig.stable
