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
sys = PSY.System(100.0, frequency = 60.0);

#Add buses
for bus in nodes_case1
    PSY.add_component!(sys, bus)
end

#Add lines
for lines in branch_case1
    PSY.add_component!(sys, lines)
end

#Add loads
for loads in loads_case1
    PSY.add_component!(sys, loads)
end

#Add infinite source
PSY.add_component!(sys, inf_gen_case1)

#Add generator
PSY.add_component!(sys, case1_gen)

##################################################
############### SOLVE PROBLEM ####################
##################################################

#Compute Y_bus after fault
sys2 = PSY.System(100.0, frequency = 60.0);
#Add buses
for bus in nodes_case1
    PSY.add_component!(sys2, bus)
end

#Add lines
for lines in branch_case1_fault
    PSY.add_component!(sys2, lines)
end
Ybus_fault = PSY.Ybus(sys2)[:, :]

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

#Solve problem in equilibrium
run_simulation!(sim, IDA(), dtmax = 0.02);

#Obtain data for angles
series = get_state_series(sim, ("Case1Gen", :δ));
series2 = get_voltagemag_series(sim, 2)
LITS.print_init_states(sim)

@test sim.solution.retcode == :Success
