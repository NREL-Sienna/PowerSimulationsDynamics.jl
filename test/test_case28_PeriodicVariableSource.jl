"""
Case 28:
This case studies a OMIB with a simple Periodic Variable Source with a single frequency component.
There is no fault, however, the system is not in steady state due to the changing source.
"""

##################################################
############### LOAD DATA ########################
##################################################

include(joinpath(TEST_FILES_DIR, "data_tests/test28.jl"))

##################################################
############### SOLVE PROBLEM ####################
##################################################

tspan = (0.0, 1.0);
step = 1e-1
tsteps = tspan[1]:step:tspan[2]

@testset "Test 28 Periodic Variable Source ResidualModel" begin
    path = (joinpath(pwd(), "test-28"))
    !isdir(path) && mkdir(path)
    try
        #Define Simulation Problem
        sim = Simulation!(
            ResidualModel,
            sys, # system
            path,
        )

        #Solve problem
        @test execute!(sim, IDA(), tspan, saveat = tsteps) == PSID.SIMULATION_FINALIZED
        results = read_results(sim)
        pvs = collect(get_components(PeriodicVariableSource, sys))[1]

        #get time domain results from PVS data
        Vt = zeros(length(tsteps)) .+ get_internal_voltage_bias(pvs)
        ω_V = get_internal_voltage_frequencies(pvs)
        for (ix, A) in enumerate(get_internal_voltage_coefficients(pvs))
            Vt += A[1] * sin.(ω_V[ix] .* tsteps) + A[2] * cos.(ω_V[ix] .* tsteps)
        end
        θt = zeros(length(tsteps)) .+ get_internal_angle_bias(pvs)
        ω_θ = get_internal_angle_frequencies(pvs)
        for (ix, A) in enumerate(get_internal_angle_coefficients(pvs))
            θt += A[1] * sin.(ω_θ[ix] .* tsteps) + A[2] * cos.(ω_θ[ix] .* tsteps)
        end

        #Obtain simulation data
        Vt_sim = get_state_series(results, ("InfBus", :Vt))
        θt_sim = get_state_series(results, ("InfBus", :θt))

        @test LinearAlgebra.norm(Vt .- Vt_sim[2]) <= 5e-3
        @test LinearAlgebra.norm(θt .- θt_sim[2]) <= 5e-3
    finally
        @info("removing test files")
        rm(path, force = true, recursive = true)
    end
end

@testset "Test 28 Periodic Variable Source MassMatrixModel" begin
    path = (joinpath(pwd(), "test-28"))
    !isdir(path) && mkdir(path)
    try
        #Define Simulation Problem
        sim = Simulation!(
            MassMatrixModel,
            sys, # system
            path,
        )

        #Solve problem
        @test execute!(sim, Rodas4(), tspan, saveat = tsteps) == PSID.SIMULATION_FINALIZED
        results = read_results(sim)
        pvs = collect(get_components(PeriodicVariableSource, sys))[1]

        #get time domain results from PVS data
        Vt = zeros(length(tsteps)) .+ get_internal_voltage_bias(pvs)
        ω_V = get_internal_voltage_frequencies(pvs)
        for (ix, A) in enumerate(get_internal_voltage_coefficients(pvs))
            Vt += A[1] * sin.(ω_V[ix] .* tsteps) + A[2] * cos.(ω_V[ix] .* tsteps)
        end
        θt = zeros(length(tsteps)) .+ get_internal_angle_bias(pvs)
        ω_θ = get_internal_angle_frequencies(pvs)
        for (ix, A) in enumerate(get_internal_angle_coefficients(pvs))
            θt += A[1] * sin.(ω_θ[ix] .* tsteps) + A[2] * cos.(ω_θ[ix] .* tsteps)
        end

        #Obtain simulation data
        Vt_sim = get_state_series(results, ("InfBus", :Vt))
        θt_sim = get_state_series(results, ("InfBus", :θt))

        @test LinearAlgebra.norm(Vt .- Vt_sim[2]) <= 5e-3
        @test LinearAlgebra.norm(θt .- θt_sim[2]) <= 5e-3
    finally
        @info("removing test files")
        rm(path, force = true, recursive = true)
    end
end
