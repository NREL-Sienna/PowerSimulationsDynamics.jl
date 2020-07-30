"""
Validation PSSE/GENROU:
This case study defines a three bus system with an infinite bus, GENROU and a load.
The fault drop the line connecting the infinite bus and GENROU.
"""

##################################################
############### LOAD DATA ########################
##################################################

include(joinpath(dirname(@__FILE__), "data_tests/test_genrou.jl"))

##################################################
############### SOLVE PROBLEM ####################
##################################################
#Define Fault: Change of YBus
Ybus_change = ThreePhaseFault(
    1.0, #change at t = 1.0
    Ybus_fault, #New YBus
) 


path = (joinpath(pwd(), "test-psse-genrou"))
!isdir(path) && mkdir(path)
try
    #Define Simulation Problem
    sim = Simulation(
        path,
        sys, #system
        (0.0, 20.0), #time span
        Ybus_change,
    ) #Type of Fault

    #Obtain small signal results for initial conditions
    small_sig = small_signal_analysis(sim)

    #Solve problem in equilibrium
    run_simulation!(sim, IDA(), dtmax = 0.001, saveat = 0.001)

    #Obtain data for angles
    series = get_state_series(sim, ("generator-102-1", :δ))
    t = series[1]
    δ = series[2]
    #Clean Extra Point at t = 1.0 from Callback
    clean_extra_timestep!(t, δ)

    series2 = get_voltagemag_series(sim, 102)

    psse_csv = joinpath(dirname(@__FILE__), "benchmarks/psse/GENROU/Test_GENROU.csv")
    t_psse, δ_psse = get_csv_delta(psse_csv)

    diff = [0.0]
    res = get_init_values_for_comparison(sim)
    for (k, v) in test_psse_genrou_init
        diff[1] += LinearAlgebra.norm(res[k] - v)
    end
    #Test Initial Condition
    @test (diff[1] < 1e-3)
    #Test Solution DiffEq
    @test sim.solution.retcode == :Success
    #Test Small Signal
    @test small_sig.stable
    #Test Transient Simulation Results
    # PSSE results are in Degrees
    @test LinearAlgebra.norm(δ - (δ_psse .* pi / 180), Inf) <= 1e-1
    #@test LinearAlgebra.norm(t - t_psat) == 0.0

finally
    @info("removing test files")
    rm(path, force = true, recursive = true)
end

