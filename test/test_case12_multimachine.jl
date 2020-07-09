"""
Case 12:
This case study a three bus system with 2 machines (Classic Model - Single Shaft: 2 State model) without loads.
The machine at bus 1 is used as a reference machine, while machine at bus 2 has a simplified droop governor (TGTypeII).
The perturbation trips one of the four (out of 5) circuits of line between buses 1 and 2, multiplying by 4 its impedance.
"""

##################################################
############### LOAD DATA ########################
##################################################

include(joinpath(dirname(@__FILE__), "data_tests/test12.jl"))

##################################################
############### SOLVE PROBLEM ####################
##################################################

#Time span
tspan = (0.0, 5.0)

#Define Fault: Change of YBus
Ybus_change = ThreePhaseFault(
    1.0, #change at t = 1.0
    Ybus_fault,
) #New YBus

sim = Simulation(
    threebus_sys, #system
    tspan, #time span
    Ybus_change, #Type of Fault
)

small_sig = small_signal_analysis(sim)

#Run simulation
run_simulation!(
    sim, #simulation structure
    IDA(),#Sundials DAE Solver
    dtmax = 0.02, #keywords arguments
)

series = get_state_series(sim, ("generator-2-1", :ω));

@test LinearAlgebra.norm(sim.x0_init - test12_x0_init) < 1e-6
@test sim.solution.retcode == :Success
@test small_sig.stable
