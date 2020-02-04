"""
Case 1:
This case study defines a classical machine against an infinite bus. The fault
drop a circuit on the (double circuit) line connecting the two buses, doubling its impedance
"""

##################################################
############### LOAD DATA ########################
##################################################

############### Data Network ########################

nodes_case1 = nodes_OMIB()

branch_case1 = branches_OMIB(nodes_case1)

#Trip of a single circuit of Line 1 -> Resistance and Reactance doubled.
branch_case1_fault = branches_OMIB_fault(nodes_case1)

loads_case1 = loads_OMIB(nodes_case1)

############### Data devices ########################

inf_gen_case1 = inf_gen_105_pu(nodes_case1)

### Case 1 Generator ###

case1_gen = dyn_gen_OMIB(nodes_case1)

######################### Dynamical System ########################

#Create system with BasePower = 100 MVA and nominal frequency 60 Hz.
sys = system_no_inv(nodes_case1, branch_case1, loads_case1, [inf_gen_case1], [case1_gen])

##################################################
############### SOLVE PROBLEM ####################
##################################################

#Compute Y_bus after fault
Ybus_fault = get_admittance_matrix(nodes_case1, branch_case1_fault)

tspan = (0.0, 30.0);

#Define Fault: Change of YBus
Ybus_change = ThreePhaseFault(
    1.0, #change at t = 1.0
    Ybus_fault,
) #New YBus

#Define Simulation Problem
sim = Simulation(
    sys, #system
    tspan, #time span
    Ybus_change,
) #Type of Fault

#Obtain small signal results for initial conditions
small_sig = small_signal_analysis(sim)

#Solve problem in equilibrium
run_simulation!(sim, IDA(), dtmax = 0.02);

#Obtain data for angles
series = get_state_series(sim, ("Case1Gen", :δ));
series2 = get_voltagemag_series(sim, 2)
LITS.print_init_states(sim)

@test sim.solution.retcode == :Success
@test small_sig.stable
