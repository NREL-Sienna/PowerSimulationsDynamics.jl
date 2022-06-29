"""
Case 38:
This case study a three bus system with one AggregateDistributedGenerationA model, one load, and one infinite source. 
TODO - describe fault 
"""
##################################################
############### LOAD DATA ########################
##################################################

include(joinpath(TEST_FILES_DIR, "data_tests/test38.jl"))

##################################################
############### SOLVE PROBLEM ####################
##################################################

tspan = (0.0, 2.0)
csv_file = joinpath(TEST_FILES_DIR, "benchmarks/psse/DERA/TEST_DERA.csv")

@testset "Test 38 DERA Model" begin
    path = (joinpath(pwd(), "test-psse-dera"))
    !isdir(path) && mkdir(path)
    try
        @warn "In test"
        for g in get_components(ThermalStandard, threebus_sys)
            case_dera = dera(g)
            #case_dera = inv_case78(g)
            add_component!(threebus_sys, case_dera, g)
        end
        display(threebus_sys)
        results = solve_powerflow(threebus_sys)
        display(results["bus_results"])
        display(results["flow_results"])
        sim = Simulation(ResidualModel, threebus_sys, path, tspan)

        y_rectangular = sim.inputs.ybus_rectangular
        Y = y_rectangular[1:3, 1:3] - 1im * y_rectangular[4:6, 1:3]
        Y_fault = copy(Y)

        Y_fault[2, 2] = Y_fault[2, 2] + 1e4 #1e9 # -2*10^9*1im 
        display(Y)
        display(Y_fault)
        ns1 = NetworkSwitch(2.0, Y_fault)
        ns2 = NetworkSwitch(2.0 + 1 / 60, Y)
        sim_residual = Simulation(
            ResidualModel,
            threebus_sys,
            path,
            tspan,
            BranchTrip(1.0, Line, "BUS 1-BUS 2-i_1"),  #[ns1,ns2],
        )

        sim_mass_matrix = Simulation(
            MassMatrixModel,
            threebus_sys,
            path,
            tspan,
            BranchTrip(1.0, Line, "BUS 1-BUS 2-i_1"),  #[ns1,ns2],
        )
        show_states_initial_value(sim)
        execute!(sim_residual, IDA(), dtmax = 0.005, saveat = 0.005)
        #execute!(sim_mass_matrix, Rodas4(),  dtmax = 0.005, saveat = 0.005)
        results = read_results(sim_residual)
        v_series = get_voltage_magnitude_series(results, 102)
        θ_series = get_voltage_angle_series(results, 102)

        p1 = plot(v_series)
        p2 = plot(θ_series)
        display(plot(p1, p2))
        Vmeas_series = get_state_series(results, ("generator-102-1", :Vmeas))
        Pmeas_series = get_state_series(results, ("generator-102-1", :Pmeas))
        Q_V_series = get_state_series(results, ("generator-102-1", :Q_V))
        Iq_series = get_state_series(results, ("generator-102-1", :Iq))
        Mult_series = get_state_series(results, ("generator-102-1", :Mult))
        Fmeas_series = get_state_series(results, ("generator-102-1", :Fmeas))
        Ip_series = get_state_series(results, ("generator-102-1", :Ip))

        M = get_csv_data(csv_file)
        t_psse, θ_psse = clean_extra_timestep!(M[:, 1], M[:, 2])

        display(θ_psse[1])
        display(θ_series[2][1] * 180 / pi)
        display(θ_psse[end])
        display(maximum(Ip_series[2]))
        display(minimum(Ip_series[2]))
        display(θ_series[2][end] * 180 / pi)
        p1 = plot(v_series, label = "Vt")
        p2 = plot(Vmeas_series)
        p3 = plot(Pmeas_series, label = "Pmeas")
        p4 = plot(Q_V_series, label = "Q_V")
        p5 = plot(Iq_series, label = "Iq")
        p6 = plot(Mult_series, label = "Mult")
        p7 = plot(Fmeas_series, label = "Fmeas")
        p8 = plot(Ip_series, label = "Ip")
        p9 = plot(θ_series, label = "θ_psid")
        p10 = plot(t_psse[400:420], θ_psse[400:420], label = "θ_psse")
        display(plot(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10))
        @warn "test warn"
    finally
        @info("removing test files")
        rm(path, force = true, recursive = true)
    end
end
