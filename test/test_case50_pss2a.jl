"""
Test for PSS model : PSS2A available in PSS/e
This case study defines a three bus system with an infinite bus in 1,
a constant impedance load in bus 2, a constant impedance load in bus 3
and a GENROE in bus 3.
The GENROE machine has the ESAC1A and the PSS2A models.
The small disturbance is a change of Vref in ESAC1A.
"""
##################################################
############### LOAD DATA ########################
##################################################

raw_file = joinpath(TEST_FILES_DIR, "benchmarks/psse/PSS2A/3_BUS_PSS2A.raw")
dyr_file = joinpath(TEST_FILES_DIR, "benchmarks/psse/PSS2A/3_BUS_PSS2A.dyr")
csv_file = joinpath(TEST_FILES_DIR, "benchmarks/psse/PSS2A/results_PSSe.csv")

function get_gen_by_number(system, number)
    for gen in get_components(Generator, system)
        if get_number(get_bus(gen)) == number
            return gen
        end
    end
end

#MACHINE GENROE

genroe_machine() = RoundRotorExponential(;
    R = 0.0,
    Td0_p = 7.50,
    Td0_pp = 0.025,
    Tq0_p = 1.50,
    Tq0_pp = 0.100,
    Xd = 2.20,
    Xq = 1.70,
    Xd_p = 0.30,
    Xq_p = 0.80,
    Xd_pp = 0.24,
    Xl = 0.15,
    Se = (0.0, 0.0),
)

#Shafts

single_mass_shaft() = PSY.SingleMass(; H = 3.50, D = 0.00)

#AVR ESAC1A

esac1a_avr() = ESAC1A(;
    Tr = 2 * eps(),
    Tb = 2 * eps(),
    Tc = 2 * eps(),
    Ka = 1000.0,
    Ta = 0.128,
    Va_lim = (min = -5.50, max = 5.50),
    Te = 2.784,
    Kf = 0.004,
    Tf = 0.864,
    Kc = 2 * eps(),
    Kd = 2 * eps(),
    Ke = 1.0,
    E_sat = (0.0, 0.0),
    Se = (0.0, 0.0),
    Vr_lim = (min = -99.0, max = 99.0),
)

#Prime Mover

fixed_torque() = PSY.TGFixed(; efficiency = 1.0)

#PSS2A

pss2a_pss() = PSY.PSS2A(;
    input_code_1 = 1,
    remote_bus_control_1 = 0,
    input_code_2 = 3,
    remote_bus_control_2 = 0,
    M_rtf = 2,
    N_rtf = 4,
    Tw1 = 1.0,
    Tw2 = 1.0,
    T6 = 0.1,
    Tw3 = 1.0,
    Tw4 = 1.0,
    T7 = 0.1,
    Ks2 = 0.001,
    Ks3 = 1.0,
    T8 = 1.0,
    T9 = 0.1,
    Ks1 = 0.75,
    T1 = 0.507,
    T2 = 0.0869,
    T3 = 0.507,
    T4 = 0.0869,
    Vst_lim = (-0.2, 0.2),
)

@testset "Test 50 PSS2A ResidualModel" begin
    path = joinpath(pwd(), "test-psse-pss2a")
    !isdir(path) && mkdir(path)
    try
        # Define system
        sys = System(raw_file, dyr_file; frequency = 50.0)

        for l in get_components(PSY.StandardLoad, sys)
            transform_load_to_constant_impedance(l)
        end

        # Get gen at bus 3
        g = get_gen_by_number(sys, 3)

        # Initially the dynamic injector is a GENCLS.
        # So it needs to be removed and replaced.

        # Get the initial dynamic injector for the gen
        # at bus 3
        dynamic_injector = get_dynamic_injector(g)

        # Remove the initial dynamic injector
        remove_component!(sys, dynamic_injector)

        # Define the new dynamic injector
        dyn_gen = DynamicGenerator(;
            name = get_name(g),
            machine = genroe_machine(),
            shaft = single_mass_shaft(),
            avr = esac1a_avr(),
            prime_mover = fixed_torque(),
            ω_ref = 1.0,
            pss = pss2a_pss(),
        )

        # Attach the new dynamic injector to the gen
        # at bus 3
        add_component!(sys, dyn_gen, g)

        # Define Perturbation
        perturbation = ControlReferenceChange(1.0, dyn_gen, :V_ref, 1.05)

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
        eigs = small_sig.eigenvalues
        @test small_sig.stable

        # Solve problem
        @test execute!(sim, IDA(); dtmax = 0.005, saveat = 0.005) ==
              PSID.SIMULATION_FINALIZED
        results = read_results(sim)

        # Obtain results
        t_psid, v1_psid = get_voltage_magnitude_series(results, 1)
        _, v2_psid = get_voltage_magnitude_series(results, 2)
        _, v3_psid = get_voltage_magnitude_series(results, 3)

        # Obtain PSSE results
        M = get_csv_data(csv_file)
        t_psse = M[:, 1]
        v1_psse = M[:, 2]
        v2_psse = M[:, 3]
        v3_psse = M[:, 4]

        # Test Transient Simulation Results
        @test LinearAlgebra.norm(t_psid - round.(t_psse, digits = 3)) == 0.0
        @test LinearAlgebra.norm(v1_psid - v1_psse, Inf) <= 1e-3
        @test LinearAlgebra.norm(v2_psid - v2_psse, Inf) <= 1e-3
        @test LinearAlgebra.norm(v3_psid - v3_psse, Inf) <= 1e-3

    finally
        @info("removing test files")
        rm(path; force = true, recursive = true)
    end
end
