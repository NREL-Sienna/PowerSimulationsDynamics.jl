"""
Case 4:
This case study a three bus system with 2 machines (Marconato: 8th order model) and an infinite source.
The fault drop the connection between buses 1 and 3, eliminating the direct connection between the infinite source
and the generator located in bus 3.
"""

##################################################
############### LOAD DATA ########################
##################################################

############### Data Network ########################

nodes_case234 = nodes_3bus()

branch_case234 = branches_3lines(nodes_case234)

#Trip of Line 1.
branch_case234_fault = branches_3lines_fault(nodes_case234)

loads_case234 = loads_3bus(nodes_case234)

############### Data devices ########################

inf_gen_case234 = inf_gen_102_pu(nodes_case234)

### Case 4 Generators ###

case4_gen2 = dyn_gen2_case4(nodes_case234)

case4_gen3 = dyn_gen3_case4(nodes_case234)

######################### Dynamical System ########################

#Create system with BasePower = 100 MVA and nominal frequency 60 Hz.
sys = system_no_inv(
    nodes_case234,
    branch_case234,
    loads_case234,
    [inf_gen_case234],
    [case4_gen2, case4_gen3],
)

##################################################
############### SOLVE PROBLEM ####################
##################################################

#Compute Y_bus after fault
Ybus_fault = get_admittance_matrix(nodes_case234, branch_case234_fault)

#time span
tspan = (0.0, 20.0);

#Initial guess
x0_guess = [
    1.02,
    1.0,
    1.0,
    0.0,
    -0.01,
    -0.01,
    -0.5, #ψq
    0.8, #ψd
    1.0, #eq_p
    0.47, #ed_p
    0.95, #eq_pp
    0.8, #ed_pp
    0.6, #δ
    1.0, #ω
    2.1, #Vf
    0.28, #Vr1
    -0.39, #Vr2,
    1.0, #Vm
    -0.7, #ψq
    0.6, #ψd
    0.81, #eq_p
    0.59, #ed_p
    0.75, #eq_pp
    0.6, #ed_pp
    0.86, #δ
    1.0, #ω
    1.7, #Vf
    0.11, #Vr1
    -0.31, #Vr2,
    1.0,
] #Vm

#Define Fault: Change of YBus
Ybus_change = ThreePhaseFault(
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

#Solve problem in equilibrium
run_simulation!(sim, IDA());

#Obtain data for angles
series = get_state_series(sim, ("Case4Gen2", :δ));

@test sim.solution.retcode == :Success
