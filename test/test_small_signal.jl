#################################################
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

sim = Simulation(sys, (0.0, 30.0))

res = small_signal_analysis(sim)

@test res.stable == true
