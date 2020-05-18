"""
Case 2:
This case study a three bus system with 2 machines (One d- One q-: 4th order model) and an infinite source.
The fault drop the connection between buses 1 and 3, eliminating the direct connection between the infinite source
and the generator located in bus 3.
"""

##################################################
############### LOAD DATA ########################
##################################################

include(joinpath(dirname(@__FILE__), "data_tests/test02.jl"))

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
    1.0, #eq_p
    0.47, #ed_p
    0.6, #δ
    1.0, #ω
    2.1, #Vf
    0.28, #Vr1
    -0.39, #Vr2,
    1.0, #Vm
    0.81, #eq_p
    0.59, #ed_p
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
    (0.0, 30.0), #time span
    Ybus_change, #Type of Fault
    initial_guess = x0_guess,
) #initial guess

small_sig = small_signal_analysis(sim)

#Solve problem in equilibrium
run_simulation!(sim, IDA())

#Obtain data for angles
series = get_state_series(sim, ("generator-2-1", :δ))

@test norm(sim.x0_init - test02_x0_init) < 1e-5
@test sim.solution.retcode == :Success
@test small_sig.stable
