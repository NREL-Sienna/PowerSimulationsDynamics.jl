using PowerSimulationsDynamics
using Sundials

"""
Case 24:
This case study a three-bus system, with two Marconato generators (in buses 1 and 2), and a constant impedance load in bus 3.
The perturbation increase the reference of mechanical power of generator-2 from 0.8 to 0.9 at t=1.0s.
"""

##################################################
############### LOAD DATA ########################
##################################################

include(joinpath(dirname(@__FILE__), "data_tests/test24.jl"))

##################################################
############### SOLVE PROBLEM ####################
##################################################

#time span
tspan = (0.0, 40.0);

#Define Fault using Callbacks
gen2 = get_dynamic_injector(get_component(Generator, sys, "generator-102-1"));
Pref_change = ControlReferenceChange(1.0, gen2, PSID.P_ref_index, 0.9);

path = (joinpath(pwd(), "test-23"))
!isdir(path) && mkdir(path)
try
    #Define Simulation Problem
    sim = Simulation!(path, sys, tspan, Pref_change)

    small_sig = small_signal_analysis(sim)

    #Solve problem in equilibrium
    execute!(sim, Sundials.IDA(), dtmax = 0.001)

    #Obtain data for angles
    series = get_state_series(sim, ("generator-102-1", :ω))

    diff = [0.0]
    res = get_init_values_for_comparison(sim)
    for (k, v) in test24_x0_init
        diff[1] += LinearAlgebra.norm(res[k] - v)
    end
    @test (diff[1] < 1e-3)
    @test sim.solution.retcode == :Success
    @test small_sig.stable
finally
    @info("removing test files")
    rm(path, force = true, recursive = true)
end
