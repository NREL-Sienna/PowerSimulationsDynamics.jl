"""
Validation PSSE/GAST:
This case study defines a three bus system with an infinite bus, GENROU and a load.
The GENROU machine has connected an ESAC1A Excitation System and a GAST TurbineGov.
The fault drop the line connecting the infinite bus and GENROU.
"""

##################################################
############### SOLVE PROBLEM ####################
##################################################

raw_file_dir = joinpath(dirname(@__FILE__), "../benchmarks/psse/GAST/14bus.raw")
dyr_file_dir = joinpath(dirname(@__FILE__), "../benchmarks/psse/GAST/14_bus_dyn_data.dyr")

raw_file_dir = joinpath(pwd(), "test/benchmarks/psse/GAST/14bus.raw")
dyr_file = joinpath(pwd(), "test/benchmarks/psse/GAST/14_bus_dyn_data.dyr")
path = pwd()
tspan = (0.0, 10.0)
function test_gast(dyr_file, csv_file, init_cond)
    path = (joinpath(pwd(), "test-psse-gast"))
    !isdir(path) && mkdir(path)
    try
        sys = System(raw_file_dir, dyr_file)

        #Compute Ybus_fault
        sys2 = System(raw_file_dir)
        #Remove line connecting bus 1 to 2.
        remove_component!(Line, sys2, "BUS 02-BUS 04-i_4")
        Ybus_fault = Ybus(sys2).data

        #Define Fault: Change of YBus
        Ybus_change = ThreePhaseFault(
            1.0, #change at t = 1.0
            Ybus_fault, #New YBus
        )

        #Define Simulation Problem
        sim = Simulation!(
            path,
            sys, #system
            tspan, #time span
            Ybus_change,
        ) #Type of Fault

        #Solve problem in equilibrium
        execute!(sim, IDA(), dtmax = 0.005, saveat = 0.005)

        #Obtain small signal results for initial conditions. Testing the simulation reset
        #small_sig = small_signal_analysis(sim; reset_simulation = true)

        #Obtain data for angles
        series = get_state_series(sim, ("generator-102-1", :δ))
        #Obtain data for voltage magnitude at bus 102
        series2 = get_voltagemag_series(sim, 102)
        t = series[1]
        δ = series[2]
        t_v = deepcopy(series2[1])
        V = series2[2]
        #Clean Extra Point at t = 1.0 from Callback
        clean_extra_timestep!(t, δ)
        clean_extra_timestep!(t_v, V)

        M = get_csv_data(csv_file)
        t_psse, δ_psse, V_psse = M[:, 1], M[:, 2], M[:, 3]

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
        #@test small_sig.stable
        #Test Transient Simulation Results
        # PSSE results are in Degrees
        @test LinearAlgebra.norm(δ - (δ_psse .* pi / 180), Inf) <= 1e-2
        @test LinearAlgebra.norm(V - V_psse, Inf) <= 1e-2
        @test LinearAlgebra.norm(t - round.(t_psse, digits = 3)) == 0.0

    finally
        @info("removing test files")
        rm(path, force = true, recursive = true)
    end
end

@testset "AC1A Tests" begin
    for (ix, name) in enumerate(names)
        @testset "$(name)" begin
            dyr_file = dyr_files[ix]
            csv_file = csv_files[ix]
            init_cond = init_conditions[ix]
            test_ac1a(dyr_file, csv_file, init_cond)
        end
    end
end
