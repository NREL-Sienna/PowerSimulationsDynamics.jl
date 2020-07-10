"""
Case 4:
This case study a three bus system with 2 machines (Marconato: 8th order model) and an infinite source.
The fault drop the connection between buses 1 and 3, eliminating the direct connection between the infinite source
and the generator located in bus 3.
"""

##################################################
############### LOAD DATA ########################
##################################################

include(joinpath(dirname(@__FILE__), "data_tests/test04.jl"))

##################################################
############### SOLVE PROBLEM ####################
##################################################

#Define Fault: Change of YBus
Ybus_change = ThreePhaseFault(
    1.0, #change at t = 1.0
    Ybus_fault,
) #New YBus

path = (joinpath(pwd(), "test-01"))
!isdir(path) && mkdir(path)
try
    #Define Simulation Problem
    sim = Simulation(
        path,
        threebus_sys, #system
        (0.0, 20.0), #time span
        Ybus_change, #Type of Fault
    ) #initial guess

    small_sig = small_signal_analysis(sim)

    #Solve problem in equilibrium
    run_simulation!(sim, IDA(), dtmax = 0.02)

    #Obtain data for angles
    series = get_state_series(sim, ("generator-2-1", :δ))

    diff = [0.0]
    res = LITS.get_dict_init_states(sim)
    for (k, v) in test04_x0_init
        diff[1] += LinearAlgebra.norm(res[k] - v)
    end
    @test (diff[1] < 1e-6)
    @test sim.solution.retcode == :Success
    @test small_sig.stable
finally
    @info("removing test files")
    rm(path, force = true, recursive = true)
end
