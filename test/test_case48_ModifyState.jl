"""
Validation PSSE/HYGOV:
This case study defines a three bus system with an infinite bus, GENROU+SEXS and a load.
The perturbation modifies the frequency of the generator by 0.01. The system is the same as Test 26,
so no eigenvalues are tested.
"""

##################################################
############### SOLVE PROBLEM ####################
##################################################

raw_file = joinpath(TEST_FILES_DIR, "benchmarks/psse/SEXS/ThreeBusMulti.raw")
dyr_file = joinpath(TEST_FILES_DIR, "benchmarks/psse/SEXS/ThreeBus_SEXS.dyr")

@testset "Test 48 ModifyState ResidualModel" begin
    path = (joinpath(pwd(), "test-48"))
    !isdir(path) && mkdir(path)
    try
        sys = System(raw_file, dyr_file)
        for l in get_components(PSY.PowerLoad, sys)
            PSY.set_model!(l, PSY.LoadModels.ConstantImpedance)
        end

        # Define Simulation Problem
        sim = Simulation!(
            ResidualModel,
            sys, #system
            path,
            (0.0, 20.0), #time span
            ModifyState(1.0, 12, 0.01), #Type of Fault
        )

        # Solve problem
        @test execute!(sim, IDA(), dtmax = 0.005, saveat = 0.005) ==
              PSID.SIMULATION_FINALIZED
        results = read_results(sim)

        # Obtain results
        t, δ = get_state_series(results, ("generator-102-1", :δ))
        _, Vt = get_voltage_magnitude_series(results, 102)
        _, ω = get_state_series(results, ("generator-102-1", :ω))

        # TODO: Test in PSSE
        @test isa(δ, Vector{Float64})
        @test isa(Vt, Vector{Float64})
        @test isa(ω, Vector{Float64})
    finally
        @info("removing test files")
        rm(path, force = true, recursive = true)
    end
end

@testset "Test 48 ModifyState MassMatrixModel" begin
    path = (joinpath(pwd(), "test-48"))
    !isdir(path) && mkdir(path)
    try
        sys = System(raw_file, dyr_file)
        for l in get_components(PSY.PowerLoad, sys)
            PSY.set_model!(l, PSY.LoadModels.ConstantImpedance)
        end

        # Define Simulation Problem
        sim = Simulation!(
            MassMatrixModel,
            sys, #system
            path,
            (0.0, 20.0), #time span
            ModifyState(1.0, 12, 0.01), #Type of Fault
        )

        # Solve problem
        @test execute!(sim, Rodas4(), dtmax = 0.005, saveat = 0.005) ==
              PSID.SIMULATION_FINALIZED
        results = read_results(sim)

        # Obtain results
        t, δ = get_state_series(results, ("generator-102-1", :δ))
        _, Vt = get_voltage_magnitude_series(results, 102)
        _, ω = get_state_series(results, ("generator-102-1", :ω))

        # TODO: Test in PSSE
        @test isa(δ, Vector{Float64})
        @test isa(Vt, Vector{Float64})
        @test isa(ω, Vector{Float64})
    finally
        @info("removing test files")
        rm(path, force = true, recursive = true)
    end
end
