"""
Case 11:
This case study a three bus system with 1 machine (One d- One q-: 4th order model), a VSM of 19 states and an infinite source. The connection between buses 2 and 3 is modeled using dynamic lines.
The perturbation trips two of the three circuits of line between buses 1 and 2, triplicating its impedance.
"""

##################################################
############### LOAD DATA ########################
##################################################

include(joinpath(dirname(@__FILE__), "data_tests/test11.jl"))

##################################################
############### SOLVE PROBLEM ####################
##################################################

#time span
tspan = (0.0, 40.0)

#Define Fault: Change of YBus
Ybus_change = ThreePhaseFault(
    1.0, #change at t = 1.0
    Ybus_fault,
) #New YBus

path = (joinpath(pwd(), "test-11"))
!isdir(path) && mkdir(path)
try
    #Define Simulation Problem
    sim = Simulation(
        path,
        threebus_sys, #system
        tspan, #time span
        Ybus_change, #Type of Fault
    )

    #Obtain small signal results for initial conditions
    small_sig = small_signal_analysis(sim)

    #Solve problem in equilibrium
    execute!(sim, IDA())

    #Obtain data for voltages
    series = get_voltagemag_series(sim, 102)
    diff = [0.0]
    res = get_init_values_for_comparison(sim)
    for (k, v) in test10_x0_init
        diff[1] += LinearAlgebra.norm(res[k] - v)
    end
    @test (diff[1] < 1e-3)
    @test sim.solution.retcode == :Success
    @test small_sig.stable
finally
    @info("removing test files")
    rm(path, force = true, recursive = true)
end
