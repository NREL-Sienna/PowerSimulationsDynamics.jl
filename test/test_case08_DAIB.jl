using LITS
using Sundials

"""
Case 8:
This case study a 19-state virtual synchronous machine against an infinite bus located at bus 1, with VSM located at bus 2.
The perturbation increase the reference power (analogy for mechanical power) from 0.5 to 0.7.
"""

##################################################
############### LOAD DATA ########################
##################################################

include(joinpath(dirname(@__FILE__), "data_tests/test08.jl"))

##################################################
############### SOLVE PROBLEM ####################
##################################################

#time span
tspan = (0.0, 4.0);

#Define Fault using Callbacks
Pref_change = LITS.ControlReferenceChange(1.0, case_inv, LITS.P_ref_index, 0.7)

path = (joinpath(pwd(), "test-08"))
!isdir(path) && mkdir(path)
try
    #Define Simulation Problem
    sim = LITS.Simulation(path, omib_sys, tspan, Pref_change)

    small_sig = small_signal_analysis(sim)

    #Solve problem in equilibrium
    run_simulation!(sim, Sundials.IDA())

    #Obtain data for angles
    series = get_state_series(sim, ("generator-102-1", :ω_oc))

    print_device_states(sim)

    diff = [0.0]
    res = LITS.get_dict_init_states(sim)
    for (k, v) in test08_x0_init
        diff[1] += LinearAlgebra.norm(res[k] - v)
    end
    @test (diff[1] < 1e-3)
    @test sim.solution.retcode == :Success
    @test small_sig.stable
finally
    @info("removing test files")
    rm(path, force = true, recursive = true)
end
