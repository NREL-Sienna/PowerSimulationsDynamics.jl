"""
Case 5:
This case study a three bus system with 1 machine located at bus 2.
The generator uses the model of a one d- one q- machine, and has a 5-mass shaft and a turbine governor.
The fault disconnects a circuit between buses 1 and 2, doubling its impedance.
"""

##################################################
############### LOAD DATA ########################
##################################################

############### Data Network ########################

nodes_case5 = nodes_3bus_case5()

branch_case5 = branches_3lines_case5(nodes_case5)

#Trip of a single circuit of Line 1 -> Resistance and Reactance doubled.
branch_case5_fault = branches_3lines_case5_fault(nodes_case5)

loads_case5 = loads_3bus_case5(nodes_case5)

############### Data devices ########################

inf_gen_case5 = inf_gen_1_pu(nodes_case5)

### Case 5 Generator ###

case5_gen = dyn_gen_case5(nodes_case5)

######################### Dynamical System ########################

#Create system with BasePower = 100 MVA and nominal frequency 60 Hz.
sys = system_no_inv(nodes_case5, branch_case5, loads_case5, [inf_gen_case5], [case5_gen])

##################################################
############### SOLVE PROBLEM ####################
##################################################

#Compute Y_bus after fault
Ybus_fault = get_admittance_matrix(nodes_case5, branch_case5_fault)

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
    1.0, #eq_p
    0.0, #ed_p
    0.05, #δ
    1.0, #ω
    0.05, #δ_HP
    1.0, #ω_HP
    0.05, #δ_IP
    1.0, #ω_LP
    0.05, #δ_LP
    1.0, #ω_LP
    0.05, #δ_Ex
    1.0, #ω_Ex
    1.0, #Vf
    0.01, #Vr1
    -0.1, #Vr2,
    1.0, #Vm
    0.0,
] #xg

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
series = get_state_series(sim, ("Case5Gen", :δ));

@test sim.solution.retcode == :Success
