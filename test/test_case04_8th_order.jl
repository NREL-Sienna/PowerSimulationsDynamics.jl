"""
Case 4:
This case study a three bus system with 2 machines (Marconato: 8th order model) and an infinite source.
The fault drop the connection between buses 1 and 3, eliminating the direct connection between the infinite source
and the generator located in bus 3.
"""

##################################################
############### LOAD DATA ########################
##################################################

include(joinpath(dirname(@__FILE__),"data_tests/test04.jl"))

##################################################
############### SOLVE PROBLEM ####################
##################################################
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
    threebus_sys, #system
    (0.0, 20.0), #time span
    Ybus_change, #Type of Fault
    initial_guess = x0_guess,
) #initial guess

#Solve problem in equilibrium
run_simulation!(sim, IDA());

#Obtain data for angles
series = get_state_series(sim, ("Case4_generator-2-1", :δ));

@test sim.solution.retcode == :Success
