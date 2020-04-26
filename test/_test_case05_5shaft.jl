
using LITS
using NLsolve
using Sundials

"""
Case 5:
This case study a three bus system with 1 machine located at bus 2.
The generator uses the model of a one d- one q- machine, and has a 5-mass shaft and a turbine governor.
The fault disconnects a circuit between buses 1 and 2, doubling its impedance.
"""


##################################################
############### LOAD DATA ########################
##################################################

include(joinpath("test","data_tests/test05.jl"))

##################################################
############### SOLVE PROBLEM ####################
##################################################

tspan = (0.0, 10.0);

#Initial guess
x0_guess = [
    1.02,
    1.0,
    1.0,
    0.0,
    -0.01,
    -0.01,
    0.0, #Gen1 δ
    1.0, #Gen1 ω
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
    threebus_sys, #system
    tspan, #time span
    Ybus_change, #Type of Fault
    initial_guess = x0_guess,
) #initial guess

#Solve problem in equilibrium
run_simulation!(sim, IDA());

#Obtain data for angles
series = get_state_series(sim, ("Case5_generator-2-1", :δ));

@test sim.solution.retcode == :Success
