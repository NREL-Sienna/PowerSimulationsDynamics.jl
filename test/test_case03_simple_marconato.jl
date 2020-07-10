"""
Case 3:
This case study a three bus system with 2 machines (Simple Marconato: 6th order model) and an infinite source.
The fault drop the connection between buses 1 and 3, eliminating the direct connection between the infinite source
and the generator located in bus 3.
"""

##################################################
############### LOAD DATA ########################
##################################################

include(joinpath(dirname(@__FILE__), "data_tests/test03.jl"))

##################################################
############### SOLVE PROBLEM ####################
##################################################

#time span
tspan = (0.0, 20.0);
#Define Fault: Change of YBus
Ybus_change = ThreePhaseFault(
    1.0, #change at t = 1.0
    Ybus_fault,
) #New YBus

#Define Simulation Problem
sim = Simulation(
    threebus_sys, #system
    tspan, #time span
    Ybus_change, #Type of Fault
)

small_sig = small_signal_analysis(sim)

#Solve problem in equilibrium
run_simulation!(sim, IDA(), dtmax = 0.02);

#Obtain data for angles
series = get_state_series(sim, ("generator-2-1", :δ));

diff = [0.0]
res = LITS.get_dict_init_states(sim)
for (k, v) in test03_x0_init
    diff[1] += LinearAlgebra.norm(res[k] - v)
end
@test (diff[1] < 1e-6)
@test sim.solution.retcode == :Success
@test small_sig.stable
