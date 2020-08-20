"""
Case GENROU + AVR Type I:
This case study defines a three bus system with an infinite bus, GENROU (with AVR TypeI) and a load.
The fault drop the line connecting the infinite bus and GENROU.
"""

##################################################
############### LOAD DATA ########################
##################################################

include(joinpath(dirname(@__FILE__), "data_tests/test17.jl"))

##################################################
############### SOLVE PROBLEM ####################
##################################################
#Define Fault: Change of YBus
Ybus_change = ThreePhaseFault(
    1.0, #change at t = 1.0
    Ybus_fault, #New YBus
)

path = (joinpath(pwd(), "test-psse-genrou-avr"))
!isdir(path) && mkdir(path)
try
    #Define Simulation Problem
    sim = Simulation(
        path,
        sys, #system
        (0.0, 30.0), #time span
        Ybus_change,
    ) #Type of Fault

    #Obtain small signal results for initial conditions
    small_sig = small_signal_analysis(sim)

    #Solve problem in equilibrium
    execute!(sim, IDA(), dtmax = 0.01)

    #Obtain data for angles
    series = get_state_series(sim, ("generator-102-1", :δ))

    diff = [0.0]
    res = get_init_values_for_comparison(sim)
    for (k, v) in test_psse_genrou_avr_init
        diff[1] += LinearAlgebra.norm(res[k] - v)
    end
    #Test Initial Condition
    @test (diff[1] < 1e-3)
    #Test Solution DiffEq
    @test sim.solution.retcode == :Success
    #Test Small Signal
    @test small_sig.stable

finally
    @info("removing test files")
    rm(path, force = true, recursive = true)
end
