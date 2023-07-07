"""
Case 51:
This case study a 16-state grid following inverter against an infinite bus located at bus 1, with the inverter located at bus 2.
The perturbation increase the reference power (analogy for mechanical power) from 0.5 to 0.7.
"""

##################################################
############### LOAD DATA ########################
##################################################

include(joinpath(TEST_FILES_DIR, "data_tests/test51.jl"))

##################################################
############### SOLVE PROBLEM ####################
##################################################

# time span
tspan = (0.0, 2.0);

#Define Fault using Callbacks
Pref_change = ControlReferenceChange(1.0, case_inv, :P_ref, 0.7)

@testset "Test 51 Grid Following Inverter ResidualModel" begin
    path = (joinpath(pwd(), "test-51"))
    !isdir(path) && mkdir(path)
    try
        # Define Simulation Problem
        sim = Simulation!(
            ResidualModel,
            omib_sys, # system
            path,
            tspan,
            Pref_change,
        )

        # Test Initial Condition
        diff_val = [0.0]
        res = get_init_values_for_comparison(sim)
        for (k, v) in test51_x0_init
            diff_val[1] += LinearAlgebra.norm(res[k] - v)
        end
        @test (diff_val[1] < 1e-3)

        # Obtain small signal results for initial conditions
        small_sig = small_signal_analysis(sim)
        eigs = small_sig.eigenvalues
        @test small_sig.stable
        # Test Eigenvalues
        @test LinearAlgebra.norm(eigs - test51_eigvals) < 1e-3

        #Solve problem in equilibrium
        @test execute!(sim, Sundials.IDA(); dtmax = 0.001, saveat = 0.005) ==
              PSID.SIMULATION_FINALIZED
        results = read_results(sim)

        power = PSID.get_activepower_series(results, "generator-102-1")
        rpower = PSID.get_reactivepower_series(results, "generator-102-1")
        ir = PSID.get_real_current_series(results, "generator-102-1")
        ii = PSID.get_imaginary_current_series(results, "generator-102-1")
        ω = PSID.get_frequency_series(results, "generator-102-1")
        @test isa(power, Tuple{Vector{Float64}, Vector{Float64}})
        @test isa(rpower, Tuple{Vector{Float64}, Vector{Float64}})
        @test isa(ir, Tuple{Vector{Float64}, Vector{Float64}})
        @test isa(ii, Tuple{Vector{Float64}, Vector{Float64}})
        @test isa(ω, Tuple{Vector{Float64}, Vector{Float64}})

    finally
        @info("removing test files")
        rm(path; force = true, recursive = true)
    end
end

@testset "Test 51 Grid Following Inverter MassMatrixModel" begin
    path = (joinpath(pwd(), "test-51"))
    !isdir(path) && mkdir(path)
    try
        # Define Simulation Problem
        sim = Simulation!(
            MassMatrixModel,
            omib_sys, # system
            path,
            tspan,
            Pref_change,
        )

        # Test Initial Condition
        diff_val = [0.0]
        res = get_init_values_for_comparison(sim)
        for (k, v) in test51_x0_init
            diff_val[1] += LinearAlgebra.norm(res[k] - v)
        end
        @test (diff_val[1] < 1e-3)

        # Obtain small signal results for initial conditions
        small_sig = small_signal_analysis(sim)
        eigs = small_sig.eigenvalues
        @test small_sig.stable
        # Test Eigenvalues
        @test LinearAlgebra.norm(eigs - test51_eigvals) < 1e-3

        #Solve problem in equilibrium
        @test execute!(sim, Rodas4(); dtmax = 0.001, saveat = 0.005) ==
              PSID.SIMULATION_FINALIZED
        results = read_results(sim)

        power = PSID.get_activepower_series(results, "generator-102-1")
        rpower = PSID.get_reactivepower_series(results, "generator-102-1")
        ir = PSID.get_real_current_series(results, "generator-102-1")
        ii = PSID.get_imaginary_current_series(results, "generator-102-1")
        ω = PSID.get_frequency_series(results, "generator-102-1")
        @test isa(power, Tuple{Vector{Float64}, Vector{Float64}})
        @test isa(rpower, Tuple{Vector{Float64}, Vector{Float64}})
        @test isa(ir, Tuple{Vector{Float64}, Vector{Float64}})
        @test isa(ii, Tuple{Vector{Float64}, Vector{Float64}})
        @test isa(ω, Tuple{Vector{Float64}, Vector{Float64}})

    finally
        @info("removing test files")
        rm(path; force = true, recursive = true)
    end
end
