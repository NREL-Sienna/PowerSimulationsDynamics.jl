"""
Case 10:
This case study a three bus system with 2 machine (Fixed EMF - Single Shaft: 2 State model), and a load in between.
The perturbation trips one of the two circuits of line between buses 2 and 3, duplicating its impedance.
"""

##################################################
############### LOAD DATA ########################
##################################################

##################################################
############### SOLVE PROBLEM ####################
##################################################

#Time span
tspan = (0.0, 30.0)

Ybus_change = ThreePhaseFault(
    1.0, #change at t = 1.0
    Ybus_fault,
); #New YBus

x0_guess = [
    1.02, #V1_R
    1.005, #V2_R
    1.0, #V3_R
    0.0, #V1_I
    0.01, #V2_I
    0.01, #V3_I
    0.1, #δ_1
    1.0, #ω_1
    0.1, #xg
    1.0, #δ_2
    0.0, #ω_2
]

sim = Simulation(
    sys, #system
    tspan, #time span
    Ybus_change, #Type of Fault
    initial_guess = x0_guess,
)

#Run simulation
run_simulation!(
    sim, #simulation structure
    IDA(),#Sundials DAE Solver
    dtmax = 0.02, #keywords arguments
)

series = get_state_series(sim, ("Case10Gen1", :ω));

@test sim.solution.retcode == :Success
