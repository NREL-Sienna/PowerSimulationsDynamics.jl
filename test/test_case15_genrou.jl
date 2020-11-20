"""
Validation PSSE/GENROU:
This case study defines a three bus system with an infinite bus, GENROU and a load.
The fault drop the line connecting the infinite bus and GENROU.
"""

##################################################
############### SOLVE PROBLEM ####################
##################################################

#Define dyr files

names = ["GENROU: Normal Saturation", "GENROU: No Saturation", "GENROU: High Saturation"]

dyr_files = [
    joinpath(dirname(@__FILE__), "benchmarks/psse/GENROU/ThreeBus_GENROU.dyr"),
    joinpath(dirname(@__FILE__), "benchmarks/psse/GENROU/ThreeBus_GENROU_NO_SAT.dyr"),
    joinpath(dirname(@__FILE__), "benchmarks/psse/GENROU/ThreeBus_GENROU_HIGH_SAT.dyr"),
]

csv_files = (
    joinpath(dirname(@__FILE__), "benchmarks/psse/GENROU/TEST_GENROU.csv"),
    joinpath(dirname(@__FILE__), "benchmarks/psse/GENROU/TEST_GENROU_NO_SAT.csv"),
    joinpath(dirname(@__FILE__), "benchmarks/psse/GENROU/TEST_GENROU_HIGH_SAT.csv"),
)

init_conditions =
    [test_psse_genrou_init, test_psse_genrou_no_sat_init, test_psse_genrou_high_sat_init]

raw_file_dir = joinpath(dirname(@__FILE__), "benchmarks/psse/GENROU/ThreeBusMulti.raw")
tspan = (0.0, 20.0)

function test_genrou(dyr_file, csv_file, init_cond)
    path = (joinpath(pwd(), "test-psse-genrou"))
    !isdir(path) && mkdir(path)
    try
        sys = System(raw_file_dir, dyr_file)

        #Define Simulation Problem
        sim = Simulation!(
            path,
            sys, #system
            tspan, #time span
            BranchTrip(1.0, "BUS 1-BUS 2-i_1"), #Type of Fault
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

        t_psse, δ_psse = get_csv_delta(csv_file)

        diff = [0.0]
        res = get_init_values_for_comparison(sim)
        for (k, v) in init_cond
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
        @test LinearAlgebra.norm(t - round.(t_psse, digits = 3)) == 0.0

    finally
        @info("removing test files")
        rm(path, force = true, recursive = true)
    end
end

@testset "GENROU Tests" begin
    for (ix, name) in enumerate(names)
        @testset "$(name)" begin
            dyr_file = dyr_files[ix]
            csv_file = csv_files[ix]
            init_cond = init_conditions[ix]
            test_genrou(dyr_file, csv_file, init_cond)
        end
    end
end
