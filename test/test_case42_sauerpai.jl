"""
Case 42:
This case study a three bus system with 2 machines (Sauer Pai: 6th order model) and an infinite source.
The fault drop the connection between buses 1 and 3, eliminating the direct connection between the infinite source
and the generator located in bus 3. 
"""

##################################################
############### LOAD DATA ########################
##################################################

include(joinpath(TEST_FILES_DIR, "data_tests/test42.jl"))

##################################################
############### SOLVE PROBLEM ####################
##################################################

# Define Fault: Change of YBus
Ybus_change = NetworkSwitch(
    1.0, #change at t = 1.0
    Ybus_fault,
) #New YBus

@testset "Test 42 SauerPai ResidualModel" begin
    path = (joinpath(pwd(), "test-42"))
    !isdir(path) && mkdir(path)
    try
        # Define Simulation Problem
        sim = Simulation!(
            ResidualModel,
            threebus_sys, #system
            path,
            (0.0, 20.0), #time span
            Ybus_change, #Type of Fault
        ) #initial guess
        dict_setpoints = get_setpoints(sim)

        for g in get_components(DynamicGenerator, threebus_sys)
            display("Power setpoints")
            display(get_P_ref(g))
        end 
        # Test Initial Condition
        diff_val = [0.0]
        res = get_init_values_for_comparison(sim)
        for (k, v) in test42_x0_init
            diff_val[1] += LinearAlgebra.norm(res[k] - v)
        end

        @test (diff_val[1] < 1e-3)

        # Obtain small signal results for initial conditions
        small_sig = small_signal_analysis(sim)
        eigs = small_sig.eigenvalues
        @test small_sig.stable

        # Test Eigenvalues
        #@test LinearAlgebra.norm(eigs - test42_eigvals) < 1e-3

        # Solve problem
        @test execute!(sim, IDA(), dtmax = 0.005, saveat = 0.005) ==
              PSID.SIMULATION_FINALIZED
        results = read_results(sim)

        # Obtain data for angles
        series = get_state_series(results, ("generator-102-1", :δ))
        t = series[1]
        δ = series[2]
        display(plot(t,δ))
    
        # Should return zeros and a warning
        series3 = get_field_current_series(results, "generator-102-1")
 
        # Obtain PSAT benchmark data
        psat_csv = joinpath(TEST_FILES_DIR, "benchmarks/psat/Test42/Test42_delta.csv")
        t_psat, δ_psat = get_csv_delta(psat_csv)

        # Test Transient Simulation Results
        @test LinearAlgebra.norm(t - t_psat) == 0.0
        @test LinearAlgebra.norm(δ - δ_psat, Inf) <= 1e-3

        power = PSID.get_activepower_series(results, "generator-102-1")
        rpower = PSID.get_reactivepower_series(results, "generator-102-1")
        @test isa(power, Tuple{Vector{Float64}, Vector{Float64}})
        @test isa(rpower, Tuple{Vector{Float64}, Vector{Float64}})
    finally
        @info("removing test files")
        rm(path, force = true, recursive = true)
    end
end

@testset "Test 42 SauerPai MassMatrixModel" begin
    path = (joinpath(pwd(), "test-04"))
    !isdir(path) && mkdir(path)
    try
        # Define Simulation Problem
        sim = Simulation!(
            MassMatrixModel,
            threebus_sys, #system,
            path,
            (0.0, 20.0), #time span
            Ybus_change, #Type of Fault
        ) #initial guess

        # Test Initial Condition
        diff_val = [0.0]
        res = get_init_values_for_comparison(sim)
        for (k, v) in test42_x0_init
            diff_val[1] += LinearAlgebra.norm(res[k] - v)
        end

        @test (diff_val[1] < 1e-3)

        # Obtain small signal results for initial conditions
        small_sig = small_signal_analysis(sim)
        eigs = small_sig.eigenvalues
        @test small_sig.stable

        # Test Eigenvalues
        @test LinearAlgebra.norm(eigs - test42_eigvals) < 1e-3

        # Solve problem
        @test execute!(sim, Rodas4(), dtmax = 0.005, saveat = 0.005) ==
              PSID.SIMULATION_FINALIZED
        results = read_results(sim)

        # Obtain data for angles
        series = get_state_series(results, ("generator-102-1", :δ))
        t = series[1]
        δ = series[2]

        # Should return zeros and a warning
        series3 = get_field_current_series(results, "generator-102-1")
        
       
        # Obtain PSAT benchmark data
        psat_csv = joinpath(TEST_FILES_DIR, "benchmarks/psat/Test42/Test42_delta.csv")
        t_psat, δ_psat = get_csv_delta(psat_csv)

        # Test Transient Simulation Results
        @test LinearAlgebra.norm(t - t_psat) == 0.0
        @test LinearAlgebra.norm(δ - δ_psat, Inf) <= 1e-3

        power = PSID.get_activepower_series(results, "generator-102-1")
        rpower = PSID.get_reactivepower_series(results, "generator-102-1")
        @test isa(power, Tuple{Vector{Float64}, Vector{Float64}})
        @test isa(rpower, Tuple{Vector{Float64}, Vector{Float64}})

    finally
        @info("removing test files")
        rm(path, force = true, recursive = true)
    end

end
