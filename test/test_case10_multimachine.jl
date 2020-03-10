"""
Case 10:
This case study a three bus system with 2 machine (Fixed EMF - Single Shaft: 2 State model), and a load in between.
The perturbation trips one of the two circuits of line between buses 2 and 3, duplicating its impedance.
"""

##################################################
############### LOAD DATA ########################
##################################################

############### Data Network ########################

nodes_case10 = nodes_multimachine()

branch_case10 = branch_multimachine(nodes_case10)

#Trip of Line 1.
branch_case10_fault = branch_multimachine_fault(nodes_case10)

loads_case10 = loads_multimachine(nodes_case10)

######## Machine Data #########

### Case 9 Generators ###
case10_gen1 = dyn_gen_case10_ref(nodes_case10)
case10_gen2 = dyn_gen_case10_2(nodes_case10)

######################### Dynamical System ########################

sys = system_no_inv_no_sources(
    nodes_case10,
    branch_case10,
    loads_case10,
    [case10_gen1, case10_gen2],
);

##################################################
############### SOLVE PROBLEM ####################
##################################################

#Compute Y_bus after fault
Ybus_fault = Ybus(branch_case10_fault, nodes_case10)[:, :]

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
