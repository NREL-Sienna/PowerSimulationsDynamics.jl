"""
Case 27:
This case study a three bus system with 1 machine (One d- One q-: 4th order model), a VSM of 19 states and an infinite source.
The test changes botht he voltage magnitude and phase angle of the source bus.
"""

##################################################
############### LOAD DATA ########################
##################################################

# Use the sme test data as Test 09
include(joinpath(dirname(@__FILE__), "data_tests/test09.jl"))

##################################################
############### SOLVE PROBLEM ####################
##################################################

####### Changing magnitude of votlage at source bus #########

#time span
tspan = (0.0, 20.0);
case_source = collect(PSY.get_components(PSY.Source, threebus_sys))[1]

#Define Fault using Callbacks
V_source_change = SourceBusVoltageChange(1.0, case_source, PSID.V_source_index, 1.1)

@testset "Test 27 Source Bus Voltage Magnitude Perturbation ImplicitModel" begin
    path = (joinpath(pwd(), "test-27"))
    !isdir(path) && mkdir(path)
    try
        #Define Simulation Problem
        sim = Simulation!(
            ImplicitModel,
            threebus_sys, # system
            path,
            tspan,
            V_source_change,
        )

        #Solve problem
        execute!(sim, IDA(), dtmax = 0.02)

        #Obtain data for angles
        series = get_state_series(sim, ("generator-103-1", :θ_oc))

        diff = [0.0]
        res = get_init_values_for_comparison(sim)
        for (k, v) in test09_x0_init
            diff[1] += LinearAlgebra.norm(res[k] - v)
        end
        @test (diff[1] < 1e-3)
        @test sim.solution.retcode == :Success
    finally
        @info("removing test files")
        rm(path, force = true, recursive = true)
    end
end

@testset "Test 27 Source Bus Voltage Magnitude Perturbation MassMatrixModel" begin
    path = (joinpath(pwd(), "test-27"))
    !isdir(path) && mkdir(path)
    try
        #Define Simulation Problem
        sim = Simulation!(
            MassMatrixModel,
            threebus_sys, # system
            path,
            tspan,
            V_source_change,
        )

        #Solve problem
        execute!(sim, Rodas5(), dtmax = 0.02)

        #Obtain data for angles
        series = get_state_series(sim, ("generator-103-1", :θ_oc))

        diff = [0.0]
        res = get_init_values_for_comparison(sim)
        for (k, v) in test09_x0_init
            diff[1] += LinearAlgebra.norm(res[k] - v)
        end
        @test (diff[1] < 1e-3)
        @test sim.solution.retcode == :Success

    finally
        @info("removing test files")
        rm(path, force = true, recursive = true)
    end
end

####### Changing angle of voltage at source bus #########

#time span
tspan = (0.0, 20.0);
case_source = collect(PSY.get_components(PSY.Source, threebus_sys))[1]

#Define Fault using Callbacks
V_source_change = SourceBusVoltageChange(1.0, case_source, PSID.θ_source_index, 0.1)

@testset "Test 27 Source Bus Voltage Angle Perturbation ImplicitModel" begin
    path = (joinpath(pwd(), "test-27"))
    !isdir(path) && mkdir(path)
    try
        #Define Simulation Problem
        sim = Simulation!(
            ImplicitModel,
            threebus_sys, # system
            path,
            tspan,
            V_source_change,
        )

        #Solve problem
        execute!(sim, IDA(), dtmax = 0.02)

        #Obtain data for angles
        series = get_state_series(sim, ("generator-103-1", :θ_oc))

        diff = [0.0]
        res = get_init_values_for_comparison(sim)
        for (k, v) in test09_x0_init
            diff[1] += LinearAlgebra.norm(res[k] - v)
        end
        @test (diff[1] < 1e-3)
        @test sim.solution.retcode == :Success
    finally
        @info("removing test files")
        rm(path, force = true, recursive = true)
    end
end

@testset "Test 27 Source Bus Voltage Angle Perturbation MassMatrixModel" begin
    path = (joinpath(pwd(), "test-27"))
    !isdir(path) && mkdir(path)
    try
        #Define Simulation Problem
        sim = Simulation!(
            MassMatrixModel,
            threebus_sys, # system
            path,
            tspan,
            V_source_change,
        )

        #Solve problem
        execute!(sim, Rodas5(), dtmax = 0.02)

        #Obtain data for angles
        series = get_state_series(sim, ("generator-103-1", :θ_oc))

        diff = [0.0]
        res = get_init_values_for_comparison(sim)
        for (k, v) in test09_x0_init
            diff[1] += LinearAlgebra.norm(res[k] - v)
        end
        @test (diff[1] < 1e-3)
        @test sim.solution.retcode == :Success

    finally
        @info("removing test files")
        rm(path, force = true, recursive = true)
    end
end
