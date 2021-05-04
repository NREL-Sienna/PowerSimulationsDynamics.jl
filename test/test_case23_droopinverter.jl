using PowerSimulationsDynamics
using Sundials

"""
Case 23:
This case study a 15-state droop grid forming inverter against an infinite bus located at bus 1, with inverter located at bus 2.
The perturbation increase the reference power (analogy for mechanical power) from 0.5 to 0.7.
"""

##################################################
############### LOAD DATA ########################
##################################################

include(joinpath(dirname(@__FILE__), "data_tests/test23.jl"))

##################################################
############### SOLVE PROBLEM ####################
##################################################

#time span
tspan = (0.0, 4.0);

#Define Fault using Callbacks
Pref_change = ControlReferenceChange(1.0, case_inv, PSID.P_ref_index, 0.7)

path = (joinpath(pwd(), "test-23"))
!isdir(path) && mkdir(path)
try
    #Define Simulation Problem
    sim = Simulation!(path, omib_sys, tspan, Pref_change)

    small_sig = small_signal_analysis(sim)

    #Solve problem in equilibrium
    execute!(sim, Sundials.IDA())

    #Obtain data for angles
    series = get_state_series(sim, ("generator-102-1", :θ_oc))

    print_device_states(sim)

    diff = [0.0]
    res = get_init_values_for_comparison(sim)
    for (k, v) in test23_x0_init
        diff[1] += LinearAlgebra.norm(res[k] - v)
    end
    @test (diff[1] < 1e-3)
    @test sim.solution.retcode == :Success
    @test small_sig.stable
finally
    @info("removing test files")
    rm(path, force = true, recursive = true)
end 