"""
Case 14:
This case study a three bus system with 1 machines (Classic Model - Single Shaft: 2 State model) and 1 Inverter without loads.
The inverter at bus 1 is used as a reference device, while machine at bus 2 has a simplified droop governor (TGTypeII).
The perturbation trips four (out of 5) circuits of line between buses 1 and 2, multiplying by 4 its impedance.
"""

##################################################
############### LOAD DATA ########################
##################################################

include(joinpath(dirname(@__FILE__), "data_tests/test14.jl"))

##################################################
############### SOLVE PROBLEM ####################
##################################################

#Time span
tspan = (0.0, 10.0)

#Define Fault: Change of YBus
Ybus_change = ThreePhaseFault(
    1.0, #change at t = 1.0
    Ybus_fault,
) #New YBus

path = (joinpath(pwd(), "test-14"))
!isdir(path) && mkdir(path)
try
    sim = Simulation(
        path,
        threebus_sys, #system
        tspan, #time span
        Ybus_change, #Type of Fault
    )

    small_sig = small_signal_analysis(sim)

    #Run simulation
    run_simulation!(
        sim, #simulation structure
        IDA(),#Sundials DAE Solver
        dtmax = 0.001, #keywords arguments
    )

    series = get_state_series(sim, ("generator-101-1", :ω_oc))

    diff = [0.0]
    res = get_init_values_for_comparison(sim)
    for (k, v) in test14_x0_init
        diff[1] += LinearAlgebra.norm(res[k] - v)
    end
    @test (diff[1] < 1e-3)
    @test sim.solution.retcode == :Success
    @test small_sig.stable
finally
    @info("removing test files")
    rm(path, force = true, recursive = true)
end
