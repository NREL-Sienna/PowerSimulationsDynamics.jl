"""
Case 9:
This case study a three bus system with 1 machine (One d- One q-: 4th order model), a VSM of 19 states and an infinite source. The connection between buses 2 and 3 is modeled using dynamic lines.
The perturbation trips two of the three circuits of line between buses 1 and 2, triplicating its impedance.
"""

##################################################
############### LOAD DATA ########################
##################################################


##################################################
############### SOLVE PROBLEM ####################
##################################################

#Initial guess
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
    1.0, #Vm
    0.5, #IL1_R
    0.5,
] #IL1_I

#Define Fault: Change of YBus
Ybus_change = LITS.ThreePhaseFault(
    1.0, #change at t = 1.0
    Ybus_fault,
) #New YBus

sim = Simulation(
    sys, #system
    tspan, #time span
    Ybus_change, #Type of perturbation
    initial_guess = x0_guess, #initial guess.
)

#Run simulation
run_simulation!(
    sim, #simulation structure
    IDA(),
) #Sundials DAE Solver

#Obtain data for voltages
series = get_voltagemag_series(sim, 2)
zoom = [
    (series[1][ix], series[2][ix])
    for (ix, s) in enumerate(series[1]) if (s > 0.90 && s < 1.6)
]

@test sim.solution.retcode == :Success
