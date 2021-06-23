"""
Case 9:
This case study a three bus system with 1 machine (One d- One q-: 4th order model), a VSM of 19 states and an infinite source.
The perturbation changes the votlage of the infinite bus to 1.1 p.u. at 1 second.
"""

##################################################
############### LOAD DATA ########################
##################################################

include(joinpath(dirname(@__FILE__), "data_tests/test09.jl"))

##################################################
############### SOLVE PROBLEM ####################
##################################################

#time span
tspan = (0.0, 20.0);
case_source=collect(PSY.get_components(PSY.Source, threebus_sys))[1]

#Define Fault using Callbacks
V_source_change = SourceBusVoltageChange(1.0, case_source, 1.1)

@testset "Test 27 Source Bus Perturbation ImplicitModel" begin
    path = (joinpath(pwd(), "test-09"))
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

@testset "Test 27 Source Bus Perturbation MassMatrixModel" begin
    path = (joinpath(pwd(), "test-09"))
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
