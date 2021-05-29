"""
Case 7:
This case study a three bus system with 2 machine located at bus 2 and 3.
The generator at bus 3 uses the model of a one d- one q- machine, and has a 5-mass shaft and a turbine governor.
The generator at bus 2 uses the model of a one d- one q- machine, and single mass shaft.
The fault disconnects a circuit between buses 1 and 2, doubling its impedance.
"""

##################################################
############### LOAD DATA ########################
##################################################

include(joinpath(dirname(@__FILE__), "data_tests/test07.jl"))

##################################################
############### SOLVE PROBLEM ####################
##################################################

tspan = (0.0, 20.0);

#Define Fault: Change of YBus
Ybus_change = NetworkSwitch(
    1.0, #change at t = 1.0
    Ybus_fault,
) #New YBus

path = (joinpath(pwd(), "test-07"))
!isdir(path) && mkdir(path)
try
    #Define Simulation Problem
    sim = Simulation!(
        ImplicitModel,
        threebus_sys, #system,
        path,
        tspan, #time span
        Ybus_change, #Type of Fault
    ) #initial guess

    small_sig = small_signal_analysis(sim)

    #Solve problem in equilibrium
    execute!(sim, IDA(), dtmax = 0.001)

    #Obtain data for angles
    series = get_state_series(sim, ("generator-103-1", :δ))
    series2 = get_state_series(sim, ("generator-103-1", :δ_hp))
    series3 = get_state_series(sim, ("generator-103-1", :δ_ip))
    series4 = get_state_series(sim, ("generator-103-1", :δ_ex))

    diff = [0.0]
    res = get_init_values_for_comparison(sim)
    for (k, v) in test07_x0_init
        diff[1] += LinearAlgebra.norm(res[k] - v)
    end
    @test (diff[1] < 1e-3)
    @test sim.solution.retcode == :Success
    @test small_sig.stable
finally
    @info("removing test files")
    rm(path, force = true, recursive = true)
end
