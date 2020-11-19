"""
Case 1:
This case study defines a classical machine against an infinite bus. The fault
drop a circuit on the (double circuit) line connecting the two buses, doubling its impedance
"""

##################################################
############### LOAD DATA ########################
##################################################

include(joinpath(dirname(@__FILE__), "data_tests/test01.jl"))

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
    sim = Simulation!(
        path,
        omib_sys, #system
        (0.0, 20.0), #time span
        Ybus_change,
    ) #Type of Fault

    #Obtain small signal results for initial conditions
    small_sig = small_signal_analysis(sim)

    #Solve problem in equilibrium
    execute!(sim, IDA(), dtmax = 0.005, saveat = 0.005)

    #Obtain data for angles
    series = get_state_series(sim, ("generator-102-1", :δ))
    t = series[1]
    δ = series[2]
    #Clean Extra Point at t = 1.0 from Callback
    clean_extra_timestep!(t, δ)

    series2 = get_voltagemag_series(sim, 102)

    #Obtain PSAT benchmark data
    psat_csv = joinpath(dirname(@__FILE__), "benchmarks/psat/Test01/Test01_delta.csv")
    psse_csv = joinpath(dirname(@__FILE__), "benchmarks/psse/Test01/Test01_delta.csv")
    t_psat, δ_psat = get_csv_delta(psat_csv)
    t_psse, δ_psse = get_csv_delta(psse_csv)

    diff = [0.0]
    res = get_init_values_for_comparison(sim)
    for (k, v) in test01_x0_init
        diff[1] += LinearAlgebra.norm(res[k] - v)
    end
    #Test Initial Condition
    @test (diff[1] < 1e-3)
    #Test Solution DiffEq
    @test sim.solution.retcode == :Success
    #Test Small Signal
    @test small_sig.stable
    #Test Transient Simulation Results
    @test LinearAlgebra.norm(t - round.(t_psse, digits = 3)) == 0.0
    # PSSE results are in Degrees
    @test LinearAlgebra.norm(δ - (δ_psse .* pi / 180), Inf) <= 2e-3
    @test LinearAlgebra.norm(t - t_psat) == 0.0
    @test LinearAlgebra.norm(δ - δ_psat, Inf) <= 1e-3

finally
    @info("removing test files")
    rm(path, force = true, recursive = true)
end
