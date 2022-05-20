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

@testset "Test 38 DERA Model" begin
    path = (joinpath(pwd(), "test-psse-dera"))
    !isdir(path) && mkdir(path)
    try
        @warn "In test"
        for g in get_components(ThermalStandard, threebus_sys)
            case_dera = dera(g)
            add_component!(threebus_sys, case_dera, g)
        end
        results = solve_powerflow(threebus_sys)
        display(results["bus_results"])
        sim = Simulation(
            ResidualModel,
            threebus_sys,
            path,
            tspan,
            BranchTrip(1.0, Line, "BUS 1-BUS 2-i_1"),
        )
        show_states_initial_value(sim)
        execute!(sim, IDA(), abstol = 1e-9)
        results = read_results(sim)
        Vmeas_series = get_state_series(results, ("generator-102-1", :Vmeas))
        Pmeas_series = get_state_series(results, ("generator-102-1", :Pmeas))
        Q_V_series = get_state_series(results, ("generator-102-1", :Q_V))
        Iq_series = get_state_series(results, ("generator-102-1", :Iq))
        Mult_series = get_state_series(results, ("generator-102-1", :Mult))
        Fmeas_series = get_state_series(results, ("generator-102-1", :Fmeas))
        Ip_series = get_state_series(results, ("generator-102-1", :Ip))

        v_series = get_voltage_magnitude_series(results, 102)
        #=         p1 = plot(v_series, label ="Vt")
                p2 = plot(Vmeas_series, label = "Vmeas")
                p3 = plot(Pmeas_series, label = "Pmeas")
                p4 = plot(Q_V_series, label = "Q_V")
                p5 = plot(Iq_series, label = "Iq")
                p6 = plot(Mult_series, label = "Mult")
                p7 = plot(Fmeas_series, label = "Fmeas")
                p8 = plot(Ip_series, label = "Ip")
                display(plot(p1, p2, p3, p4, p5, p6, p7,p8 )) =#
    finally
        @info("removing test files")
        rm(path, force = true, recursive = true)
    end
end
