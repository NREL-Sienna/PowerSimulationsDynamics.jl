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

#Define Fault: Change of YBus
Ybus_change = ThreePhaseFault(
    1.0, #change at t = 1.0
    Ybus_fault,
) #New YBus

path = (joinpath(pwd(), "test-02"))
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
    execute!(sim, IDA(), dtmax = 0.005, saveat = 0.005)

    #Obtain data for angles
    series = get_state_series(sim, ("generator-102-1", :δ))
    t = series[1]
    δ = series[2]

    #Obtain PSAT benchmark data
    psat_csv = joinpath(dirname(@__FILE__), "benchmarks/psat/Test02/Test02_delta.csv")
    t_psat, δ_psat = get_csv_delta(psat_csv)

    #Clean Extra Point at t = 1.0 from Callback
    clean_extra_timestep!(t, δ)

    diff = [0.0]
    res = get_init_values_for_comparison(sim)
    for (k, v) in test02_x0_init
        diff[1] += LinearAlgebra.norm(res[k] - v)
    end
    @test (diff[1] < 1e-3)
    @test sim.solution.retcode == :Success
    @test small_sig.stable
    @test LinearAlgebra.norm(t - t_psat) == 0.0
    @test LinearAlgebra.norm(δ - δ_psat, Inf) <= 1e-3
finally
    @info("removing test files")
    rm(path, force = true, recursive = true)
end
