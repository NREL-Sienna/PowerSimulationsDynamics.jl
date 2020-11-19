"""
Validation PSSE/TGOV1:
This case study defines a three bus system with an infinite bus, GENROU+AC1A+TGOV1 and a load.
The fault drop the line connecting the infinite bus and GENROU.
"""

##################################################
############### SOLVE PROBLEM ####################
##################################################

raw_file = joinpath(dirname(@__FILE__), "benchmarks/psse/TGOV1/ThreeBusMulti.raw")
dyr_file = joinpath(dirname(@__FILE__), "benchmarks/psse/TGOV1/ThreeBus_TGOV1.dyr")
csv_file = joinpath(dirname(@__FILE__), "benchmarks/psse/TGOV1/TEST_TGOV1.csv")

#Construct system
sys = System(raw_file, dyr_file);

#Construct fault
sys2 = System(raw_file)
line_name = "BUS 1-BUS 2-i_1"
remove_component!(Line, sys2, line_name)
Ybus_fault = Ybus(sys2).data;
Ybus_change = ThreePhaseFault(
    1.0, #change at t = 1.0
    Ybus_fault, #New YBus
);

path = (joinpath(pwd(), "test-psse-tgov1"))
!isdir(path) && mkdir(path)
try
    sim = Simulation!(
        path,
        sys, #system
        (0.0, 20.0), #time span
        Ybus_change, #Type of Fault
    )

    #Solve problem in equilibrium
    execute!(sim, IDA(), dtmax = 0.005, saveat = 0.005)

    #Obtain small signal results for initial conditions
    #NOT WORKING DUE TO TYPES ON EXECUTE
    #small_sig = small_signal_analysis(sim)

    series = get_state_series(sim, ("generator-102-1", :δ))
    t = series[1]
    δ = series[2]
    #Clean Extra Point at t = 1.0 from Callback
    clean_extra_timestep!(t, δ)

    #Obtain PSSE results
    t_psse, δ_psse = get_csv_delta(csv_file)

    diff = [0.0]
    res = get_init_values_for_comparison(sim)
    for (k, v) in test_psse_tgov1_init
        diff[1] += LinearAlgebra.norm(res[k] - v)
    end
    #Test Initial Condition
    @test (diff[1] < 1e-3)
    #Test Solution DiffEq
    @test sim.solution.retcode == :Success
    #Test Small Signal
    #@test small_sig.stable
    #Test Transient Simulation Results
    # PSSE results are in Degrees
    @test LinearAlgebra.norm(δ - (δ_psse .* pi / 180), Inf) <= 1e-2
    @test LinearAlgebra.norm(t - round.(t_psse, digits = 3)) == 0.0
finally
    @info("removing test files")
    rm(path, force = true, recursive = true)
end
