"""
Validation Eigenvalues with Andes
This case study defines an 11-bus that uses GENROU, SEXS and TGOV1.
We compare the obtained eigenvalues with Python based tool ANDES.
Results are slighty different due to the difference on shaft and
turbine governor power-torque approximations.
"""

##################################################
############### SOLVE PROBLEM ####################
##################################################

raw_file = joinpath(TEST_FILES_DIR, "benchmarks/andes/test36/11BUS_KUNDUR.raw")
dyr_file = joinpath(TEST_FILES_DIR, "benchmarks/andes/test36/11BUS_KUNDUR_TGOV.dyr")
eigs_andes_csv = joinpath(TEST_FILES_DIR, "benchmarks/andes/test36/eigs_tgov_andes.csv")

@testset "Test 36 Eigenvalues" begin
    path = (joinpath(pwd(), "test-eigs"))
    !isdir(path) && mkdir(path)
    try
        sys = System(raw_file, dyr_file)

        # Define Simulation Problem
        sim = Simulation!(
            ResidualModel,
            sys, #system
            path,
            (0.0, 20.0), #time span
        )
        ss = small_signal_analysis(sim)
        eigs = ss.eigenvalues

        A = strip.(readdlm(eigs_andes_csv, ','))
        eigs_andes = parse.(Complex{Float64}, A)
        @test LinearAlgebra.norm(eigs - eigs_andes, 2) / length(eigs) < 1.0
    finally
        @info("removing test files")
        rm(path, force = true, recursive = true)
    end
end
