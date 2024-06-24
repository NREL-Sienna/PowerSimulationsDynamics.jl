"""
Case 1:
This case study defines a classical machine against an infinite bus. Sensitivitiy
Analysis is performed for the inertia of the generator. The fault is a change in the 
source voltage. This test also checks that incorrect user data is handled correctly. 
"""
##################################################
############### LOAD DATA ######################                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               ##
##################################################

omib_sys = build_system(PSIDTestSystems, "psid_test_omib")
ieee_9bus_sys = build_system(PSIDTestSystems, "psid_test_ieee_9bus")
andes_11bus_sys = build_system(PSIDSystems, "psid_11bus_andes")
##################################################
############### SOLVE PROBLEM ####################
##################################################

s_device = get_component(Source, omib_sys, "InfBus")
s_change = SourceBusVoltageChange(1.0, s_device, :V_ref, 1.02)
using PlotlyJS

#NOTES ON SENSITIVITY ALGORITHMS FROM SCIMLSENSITIVITY               
#ReverseDiffVJP and EnzymeVJP only options compatible with Hybrid DEs (DEs with callbacks)
#BacksolveAdjoint prone to instabilities whenever the Lipschitz constant is sufficiently large (stiff equations, PDE discretizations, and many other contexts) 
#For ForwardDiffSensitivity, convert_tspan=true is needed for hybrid equations. 
@testset "Test Gradients - OMIB; H" begin
    path = mktempdir()
    try
        sim = Simulation!(
            MassMatrixModel,
            omib_sys,
            path,
            (0.0, 5.0),
            s_change,
        )

        #GET GROUND TRUTH DATA 
        execute!(sim, Rodas4(); abstol = 1e-9, reltol = 1e-9, dtmax = 0.005, saveat = 0.005)
        res = read_results(sim)
        t, δ_gt = get_state_series(res, ("generator-102-1", :δ))

        #GET PARAMETER VALUES 
        p = get_parameter_values(sim, [("generator-102-1", :Shaft, :H)])

        function f_loss(δ, δ_gt)
            #display(plot([scatter(; y = δ_gt), scatter(; y = δ)]))
            return sum(abs.(δ - δ_gt))
        end

        #GET SENSITIVITY FUNCTIONS 
        f_forward, f_grad = get_sensitivity_functions(
            sim,
            [("generator-102-1", :Shaft, :H)],
            [("generator-102-1", :δ)],
            Rodas4(),
            f_loss;
            sensealg = ForwardDiffSensitivity(),
            abstol = 1e-9,
            reltol = 1e-9,
            dtmax = 0.005,
            saveat = 0.005,
        )

        loss_zero = f_forward(p, δ_gt)
        loss_non_zero_1 = f_forward([3.2], δ_gt)
        loss_non_zero_2 = f_forward(p, δ_gt .* 2)
        @test loss_zero == 0.0
        @test loss_non_zero_1 != 0.0
        @test loss_non_zero_2 != 0.0
        grad_zero = f_grad(p, δ_gt)
        grad_nonzero_1 = f_grad([3.14], δ_gt)
        grad_nonzero_2 = f_grad([3.15], δ_gt)
        @test isapprox(grad_zero[1], 0.0, atol = 1.0)
        @test isapprox(grad_nonzero_1[1], -8.0, atol = 1.0)
        @test isapprox(grad_nonzero_2[1], 8.0; atol = 1.0)
    finally
        @info("removing test files")
        rm(path; force = true, recursive = true)
    end
end

@testset "Test Optimization - OMIB; H" begin
    path = mktempdir()
    try
        sim = Simulation!(
            MassMatrixModel,
            omib_sys,
            path,
            (0.0, 5.0),
            s_change,
        )

        #GET GROUND TRUTH DATA 
        execute!(sim, Rodas4(); abstol = 1e-9, reltol = 1e-9, dtmax = 0.005, saveat = 0.005)
        res = read_results(sim)
        t, δ_gt = get_state_series(res, ("generator-102-1", :δ))

        #GET PARAMETER VALUES 
        p = get_parameter_values(sim, [("generator-102-1", :Shaft, :H)])

        function f_loss(δ, δ_gt)
            #display(plot([scatter(; y = δ_gt), scatter(; y = δ)]))
            return sum(abs.(δ - δ_gt))
        end

        #GET SENSITIVITY FUNCTIONS 
        f_forward, f_grad = get_sensitivity_functions(
            sim,
            [("generator-102-1", :Shaft, :H)],
            [("generator-102-1", :δ)],
            Rodas4(),
            f_loss;
            sensealg = ForwardDiffSensitivity(),
            abstol = 1e-9,
            reltol = 1e-9,
            dtmax = 0.005,
            saveat = 0.005,
        )
        #H_values = [] 
        #loss_values = []
        function callback(u, l)
            #push!(H_values, u.u[1])
            #push!(loss_values, l)
            return false
        end
        optfun = OptimizationFunction{false}(
            (u, p) -> f_forward(u, δ_gt);
            grad = (res, u, p) -> res .= f_grad(u, δ_gt),
        )
        optprob = OptimizationProblem{false}(optfun, [3.14])
        sol = Optimization.solve(
            optprob,
            OptimizationOptimisers.Adam(0.002);
            callback = callback,
            maxiters = 3,
        )
        @test isapprox(sol.u[1], 3.144000187737475, atol = 1e-6)
    finally
        @info("removing test files")
        rm(path; force = true, recursive = true)
    end
end

@testset "Test Gradients - OMIB; Xd_p" begin
    path = mktempdir()
    try
        sim = Simulation!(
            MassMatrixModel,
            omib_sys,
            path,
            (0.0, 5.0),
            s_change,
        )

        #GET GROUND TRUTH DATA 
        execute!(sim, Rodas4(); abstol = 1e-9, reltol = 1e-9, dtmax = 0.005, saveat = 0.005)
        res = read_results(sim)
        t, δ_gt = get_state_series(res, ("generator-102-1", :δ))

        #GET PARAMETER VALUES 
        p = get_parameter_values(sim, [("generator-102-1", :Machine, :Xd_p)])

        function plot_traces(δ, δ_gt)
            display(plot([scatter(; y = δ_gt), scatter(; y = δ)]))
        end
        EnzymeRules.inactive(::typeof(plot_traces), args...) = nothing
        function f_loss(δ, δ_gt)
            #plot_traces(δ, δ_gt)
            return sum(abs.(δ - δ_gt))
        end

        #GET SENSITIVITY FUNCTIONS 
        f_forward, f_grad = get_sensitivity_functions(
            sim,
            [("generator-102-1", :Machine, :Xd_p)],
            [("generator-102-1", :δ)],
            Rodas4(),
            f_loss;
            sensealg = ForwardDiffSensitivity(),
            abstol = 1e-9,
            reltol = 1e-9,
            dtmax = 0.005,
            saveat = 0.005,
        )

        loss_zero = f_forward(p, δ_gt)
        loss_non_zero_1 = f_forward(p * 1.01, δ_gt)
        loss_non_zero_2 = f_forward(p * 0.99, δ_gt)
        @test isapprox(loss_zero, 0.0, atol = 1e-9)
        @test loss_non_zero_1 != 0.0
        @test loss_non_zero_2 != 0.0
        grad_zero = f_grad(p, δ_gt)
        grad_nonzero_1 = f_grad(p * 1.01, δ_gt)
        grad_nonzero_2 = f_grad(p * 0.99, δ_gt)
        @test isapprox(grad_zero[1], 497, atol = 1.0)
        @test isapprox(grad_nonzero_1[1], 499.0, atol = 1.0)
        @test isapprox(grad_nonzero_2[1], -498.0; atol = 1.0)
    finally
        @info("removing test files")
        rm(path; force = true, recursive = true)
    end
end
