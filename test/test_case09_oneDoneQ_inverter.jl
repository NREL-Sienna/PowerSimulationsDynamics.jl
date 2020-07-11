"""
Case 9:
This case study a three bus system with 1 machine (One d- One q-: 4th order model), a VSM of 19 states and an infinite source.
The perturbation increase the reference power (analogy for mechanical power) of the inverter from 1.0 to 1.2.
"""

##################################################
############### LOAD DATA ########################
##################################################

include(joinpath(dirname(@__FILE__), "data_tests/test09.jl"))

##################################################
############### SOLVE PROBLEM ####################
##################################################

#time span
tspan = (0.0, 20.0);
case_inv = collect(PSY.get_components(PSY.DynamicInverter, threebus_sys))[1]

#Define Fault using Callbacks
Pref_change = LITS.ControlReferenceChange(1.0, case_inv, LITS.P_ref_index, 1.2)

path = (joinpath(pwd(), "test-09"))
!isdir(path) && mkdir(path)
try
    #Define Simulation Problem
    sim = LITS.Simulation(path, threebus_sys, tspan, Pref_change)

    small_sig = small_signal_analysis(sim)

    #Solve problem in equilibrium
    run_simulation!(sim, IDA(), dtmax = 0.02)

    #Obtain data for angles
    series = get_state_series(sim, ("generator-103-1", :θ_oc))

    diff = [0.0]
    res = LITS.get_dict_init_states(sim)
    for (k, v) in test09_x0_init
        diff[1] += LinearAlgebra.norm(res[k] - v)
    end
    @test (diff[1] < 1e-6)
    @test sim.solution.retcode == :Success
    @test small_sig.stable
finally
    @info("removing test files")
    rm(path, force = true, recursive = true)
end
