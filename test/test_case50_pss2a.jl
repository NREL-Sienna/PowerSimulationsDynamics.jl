"""
Test for PSS model : PSS2A available in PSS/e
This case study defines a two bus system with a GENROU machine in bus 1
and an infinite bus in 2.
The GENROU machine has the SEXS and the PSS2A models.
The small disturbance is a change of Vref in SEXS.
"""
##################################################
############### LOAD DATA ########################
##################################################

raw_file = joinpath(TEST_FILES_DIR, "benchmarks/psse/PSS2A/OMIB.raw")
dyr_file = joinpath(TEST_FILES_DIR, "benchmarks/psse/PSS2A/OMIB_GENCLS.dyr")
csv_file = joinpath(TEST_FILES_DIR, "benchmarks/psse/PSS2A/results_PSSe.csv")

function get_gen_by_number(system, number)
    for gen in get_components(Generator, system)
        if get_number(get_bus(gen)) == number
            return gen
        end
    end
end

#MACHINE GENROU

genrou_machine() = PSY.RoundRotorQuadratic(;
    R = 0.0,
    Td0_p = 7.0,
    Td0_pp = 999.0,
    Tq0_p = 0.4,
    Tq0_pp = 999.0,
    Xd = 2.20,
    Xq = 2.20,
    Xd_p = 0.30,
    Xq_p = 0.30,
    Xd_pp = 0.30,
    Xl = 0.2,
    Se = (0.0, 0.0),
)

#Shafts

single_mass_shaft() = PSY.SingleMass(H = 4.00, D = 0.00)

#AVR SEXS

sexs_avr() = PSY.SEXS(;
    Ta_Tb = 1.0,
    Tb = 1.0,
    K = 120.0,
    Te = 0.4,
    V_lim = (min = -10.0, max = 10.0),
)

#Prime Mover

fixed_torque() = PSY.TGFixed(efficiency = 1.0)

#PSS2A

pss2a_pss() = PSY.PSS2A(;
    input_code_1 = 1,
    remote_bus_control_1 = 0,
    input_code_2 = 3,
    remote_bus_control_2 = 0,
    M_rtf = 5,
    N_rtf = 1,
    Tw1 = 1.5,
    Tw2 = 1.5,
    T6 = 0.0,
    Tw3 = 1.5,
    Tw4 = 0.0,
    T7 = 1.5,
    Ks2 = 0.1875,
    Ks3 = 1.0,
    T8 = 0.5,
    T9 = 0.1,
    Ks1 = 2.0,
    T1 = 0.59451,
    T2 = 0.0447,
    T3 = 0.59451,
    T4 = 0.0447,
    Vst_lim = (-0.1, 0.1),
)

@testset "Test 50 PSS2A ResidualModel" begin
    path = joinpath(pwd(), "test-psse-pss2a")
    !isdir(path) && mkdir(path)
    try
        # Define system
        sys = System(raw_file, dyr_file)

        for l in get_components(PSY.StandardLoad, sys)
            transform_load_to_constant_impedance(l)
        end

        # Get gen at bus 1
        g = get_gen_by_number(sys, 1)

        # Initially the dynamic injector is a GENCLS.
        # So it needs to be removed and replaced.

        # Get the initial dynamic injector for the gen
        # at bus 1
        dynamic_injector = get_dynamic_injector(g)

        # Remove the initial dynamic injector
        remove_component!(sys, dynamic_injector)

        # Define the new dynamic injector
        dyn_gen = DynamicGenerator(;
            name = get_name(g),
            machine = genrou_machine(),
            shaft = single_mass_shaft(),
            avr = sexs_avr(),
            prime_mover = fixed_torque(),
            ω_ref = 1.0,
            pss = pss2a_pss(),
        )

        # Attach the new dynamic injector to the gen
        # at bus 1
        add_component!(sys, dyn_gen, g)

        # Define Perturbation
        perturbation = ControlReferenceChange(1.0, dyn_gen, :V_ref, 1.04691)

        # Define Simulation Problem
        sim = Simulation(
            ResidualModel,
            sys, #system
            path,
            (0.0, 20.0), #time span
            perturbation, #Type of Fault
        )

        # Test Initial Condition
        diff_val = [0.0]
        res = get_init_values_for_comparison(sim)
        for (k, v) in test50_x0_init
            diff_val[1] += LinearAlgebra.norm(res[k] - v)
        end

        @test diff_val[1] < 1e-3

        # Obtain small signal results for initial conditions
        small_sig = small_signal_analysis(sim)
        # @test small_sig.stable

        # In this test several blocks are bypassed resulting in problems for
        # the small signal stability analysis at initial conditions.

        # Solve problem
        @test execute!(sim, IDA(); dtmax = 0.005, saveat = 0.005) ==
              PSID.SIMULATION_FINALIZED
        results = read_results(sim)

        # Obtain results
        t_psid, Pe_psid = get_activepower_series(results, "generator-1-1")
        _, v1_psid = get_voltage_magnitude_series(results, 1)
        _, omega_psid = get_state_series(results, ("generator-1-1", :ω))

        # Obtain PSSE results
        M = get_csv_data(csv_file)
        t_psse = M[:, 1]
        Pe_psse = M[:, 2]
        v1_psse = M[:, 3]
        omega_psse = M[:, 4]

        # Test Transient Simulation Results
        @test LinearAlgebra.norm(t_psid - round.(t_psse, digits = 3)) == 0.0
        @test LinearAlgebra.norm(Pe_psid - Pe_psse, Inf) <= 5e-3
        @test LinearAlgebra.norm(v1_psid - v1_psse, Inf) <= 2e-3
        @test LinearAlgebra.norm(omega_psid - omega_psse, Inf) <= 1e-4

    finally
        @info("removing test files")
        rm(path; force = true, recursive = true)
    end
end
