"""
Test for PSS model : STAB1 available in PSS/e
This case study defines a two bus system with an infinite bus in 2,
and a GENSAL in bus 1.
The GENSAL machine has the SEXS and the STAB1 models
The small disturbance is a change of Vref in SEXS
The system is small-signal stable thanks to the PSS
"""
##################################################
############### LOAD DATA ########################
##################################################

raw_file = joinpath(TEST_FILES_DIR, "benchmarks/psse/STAB1/OMIB_SSS.raw")
dyr_file = joinpath(TEST_FILES_DIR, "benchmarks/psse/STAB1/OMIB_SSS.dyr")

sys = System(raw_file, dyr_file)
for l in get_components(PSY.PowerLoad, sys)
    PSY.set_model!(l, PSY.LoadModels.ConstantImpedance)
end

##################################################
############### SOLVE PROBLEM ####################
##################################################

csv_file = joinpath(TEST_FILES_DIR, "benchmarks/psse//STAB1/results_PSSe.csv")

@testset "Test 41 STAB1 ResidualModel" begin
    path = joinpath(pwd(), "test-psse-stab1")
    !isdir(path) && mkdir(path)
    try
        # Define Simulation Problem

        gen = first(get_components(Generator, sys))
        dynamic_injector = get_dynamic_injector(gen)

        for g in get_components(Generator, sys)
            #Find the generator at bus 1
            if get_number(get_bus(g)) == 1
                gen = g
                dynamic_injector = get_dynamic_injector(g)
            end
        end

        perturbation = ControlReferenceChange(1.0, dynamic_injector, :V_ref, 1.0472)

        sim = Simulation(
            ResidualModel,
            sys,
            path,
            (0.0, 20.0), #time span
            perturbation,
        )

        # Test Initial Condition
        diff_val = [0.0]
        res = get_init_values_for_comparison(sim)
        for (k, v) in test41_x0_init
            diff_val[1] += LinearAlgebra.norm(res[k] - v)
        end

        @test diff_val[1] < 1e-3

        # Obtain small signal results for initial conditions
        small_sig = small_signal_analysis(sim)
        eigs = small_sig.eigenvalues
        @test small_sig.stable

        # Test Eigenvalues
        @test LinearAlgebra.norm(eigs - test41_eigvals) < 1e-3

        # Solve problem
        @test execute!(sim, IDA(), dtmax = 0.005, saveat = 0.005) ==
              PSID.SIMULATION_FINALIZED
        results = read_results(sim)

        # Obtain results

        t_psid, v1_psid = get_voltage_magnitude_series(results, 1)
        _, ω_psid = get_state_series(results, ("generator-1-1", :ω))

        # Obtain PSSE results
        M = get_csv_data(csv_file)

        t_psse = M[:, 1]
        v1_psse = M[:, 2]
        ω_psse = M[:, 3] .+ 1.0

        # Test Transient Simulation Results

        @test LinearAlgebra.norm(t_psid - round.(t_psse, digits = 3)) == 0.0
        @test LinearAlgebra.norm(v1_psid - v1_psse, Inf) <= 1e-3
        @test LinearAlgebra.norm(ω_psid - ω_psse, Inf) <= 1e-3

    finally
        @info("removing test files")
        rm(path, force = true, recursive = true)
    end
end
