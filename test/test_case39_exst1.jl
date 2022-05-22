"""
Test for AVR model : EXST1 available in PSS/e
This case study defines a four bus system with an infinite bus in 1, 
a GENSAL in bus 2 and a constant impedance load in bus 3
The GENSAL machine has the EXST1 and a HYGOV.
The disturbance is the outage of one line between buses 1 and 4
"""

##################################################
############### SOLVE PROBLEM ####################
##################################################

raw_file = joinpath(TEST_FILES_DIR, "benchmarks/psse/EXST1/TVC_System_32.raw")
dyr_file = joinpath(TEST_FILES_DIR, "benchmarks/psse/EXST1/TVC_System.dyr")
csv_file = joinpath(TEST_FILES_DIR, "benchmarks/psse//EXST1/results_PSSe.csv")


@testset "Test 39 EXST1 ResidualModel" begin
    path = (joinpath(pwd(), "test-psse-exst1"))
    !isdir(path) && mkdir(path)
    try
        sys = System(raw_file, dyr_file)
        source = first(get_components(Source, sys))
        PSY.set_X_th!(source, 0.01)

        for l in get_components(PSY.PowerLoad, sys)
            PSY.set_model!(l, PSY.LoadModels.ConstantImpedance)
        end

        gen = first(get_components(Generator, sys))
        dynamic_injector = get_dynamic_injector(gen)

        for g in get_components(Generator, sys)
            #Find the generator at bus 2
            if get_number(get_bus(g)) == 2
                gen = g
                dynamic_injector = get_dynamic_injector(g)
            end
        end

        machine = deepcopy(get_machine(dynamic_injector))
        shaft = deepcopy(get_shaft(dynamic_injector))
        tg = deepcopy(get_prime_mover(dynamic_injector))
        pss = deepcopy(get_pss(dynamic_injector))
        avr = deepcopy(get_avr(dynamic_injector))

        new_dynamic_injector = DynamicGenerator(name = get_name(gen), machine = machine, shaft=shaft, avr = avr_exst1(), prime_mover = tg, ω_ref = 1.0, pss = pss)


        # Define Simulation Problem
        sim = Simulation(
            ResidualModel,
            sys, #system
            path,
            (0.0, 20.0), #time span
            BranchTrip(1.0, Line, "BUS1-BUS4-i_1"), #Type of Fault
        )
        
    
        # Test Initial Condition
        diff = [0.0]
        res = get_init_values_for_comparison(sim)
        for (k, v) in test39_x0_init
            diff[1] += LinearAlgebra.norm(res[k] - v)
        end

        @test (diff[1] < 1e-3)

        # Obtain small signal results for initial conditions
        small_sig = small_signal_analysis(sim)
        eigs = small_sig.eigenvalues
        @test small_sig.stable

        # Test Eigenvalues
        @test LinearAlgebra.norm(eigs - test39_eigvals) < 1e-3

        # Solve problem
        @test execute!(sim, IDA(), dtmax = 0.005, saveat = 0.005) ==
              PSID.SIMULATION_FINALIZED
        results = read_results(sim)

        # Obtain results

        t_psid, v2_psid = get_voltage_magnitude_series(results, 2)
        _, v3_psid = get_voltage_magnitude_series(results, 3)
        _, ω_psid = get_state_series(results, ("generator-2-1", :ω));

         
        # Obtain PSSE results
        M = get_csv_data(csv_file)
        t_psse = M[:, 1] 
        v1_psse = M[:, 2] 
        v2_psse = M[:, 3]
        v3_psse = M[:, 4] 
        v4_psse = M[:, 5] 
        ω_psse = M[:, 6] .+ 1.0
        efd_psse = M[:, 7]

      

        # Test Transient Simulation Results

        @test LinearAlgebra.norm(t_psid - round.(t_psse, digits = 3)) == 0.0
        @test LinearAlgebra.norm(v2_psid - v2_psse, Inf) <= 1e-3
        @test LinearAlgebra.norm(v3_psid - v3_psse, Inf) <= 1e-3
        @test LinearAlgebra.norm(ω_psid - ω_psse, Inf) <= 1e-3

   
    finally
        @info("removing test files")
        rm(path, force = true, recursive = true)
    end
end


