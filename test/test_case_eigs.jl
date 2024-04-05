"""
Validation Eigenvalues with Andes
This case study defines an 11-bus that uses GENROU, SEXS and TGOV1.
We compare the obtained eigenvalues with Python based tool ANDES.
Results are slighty different due to the difference on shaft and
turbine governor power-torque approximations. Removing the division
for ω for output torque in TGOV1 model produces exact eigenvalue results.
"""

##################################################
############### SOLVE PROBLEM ####################
##################################################

eigs_andes_csv = joinpath(TEST_FILES_DIR, "benchmarks/andes/test36/eigs_tgov_andes.csv")

@testset "Test 36 Eigenvalues" begin
    path = mktempdir()
    try
        sys = build_system(PSIDSystems, "psid_11bus_andes")
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
        CRC.@ignore_derivatives @info("removing test files")
        rm(path; force = true, recursive = true)
    end
end
