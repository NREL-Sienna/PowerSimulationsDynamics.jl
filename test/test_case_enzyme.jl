"""
Case 1:
This case study defines a classical machine against an infinite bus. Sensitivitiy
Analysis is performed for the inertia of the generator. The fault is a change in the 
source voltage. This test also checks that incorrect user data is handled correctly. 
"""
##################################################
############### LOAD DATA ######################                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               ##
##################################################

ieee_9bus_sys = build_system(PSIDTestSystems, "psid_test_ieee_9bus")
andes_11bus_sys = build_system(PSIDSystems, "psid_11bus_andes")
##################################################
############### SOLVE PROBLEM ####################
##################################################

using PlotlyJS

#NOTES ON SENSITIVITY ALGORITHMS FROM SCIMLSENSITIVITY               
#ReverseDiffVJP and EnzymeVJP only options compatible with Hybrid DEs (DEs with callbacks)
#BacksolveAdjoint prone to instabilities whenever the Lipschitz constant is sufficiently large (stiff equations, PDE discretizations, and many other contexts) 
#For ForwardDiffSensitivity, convert_tspan=true is needed for hybrid equations. 
@testset "Test Gradients - OMIB; H; SourceBusVoltageChange" begin
    path = mktempdir()
    try
        omib_sys = build_system(PSIDTestSystems, "psid_test_omib")
        s_device = get_component(Source, omib_sys, "InfBus")
        s_change = SourceBusVoltageChange(1.0, s_device, :V_ref, 1.02)
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

        function plot_traces(δ, δ_gt)
            display(plot([scatter(; y = δ_gt), scatter(; y = δ)]))
        end
        EnzymeRules.inactive(::typeof(plot_traces), args...) = nothing
        function f_loss(states, δ_gt)
            #plot_traces(states[1], δ_gt)
            return sum(abs.(states[1] - δ_gt))
        end

        #GET SENSITIVITY FUNCTIONS 
        f_forward, f_grad, _ = get_sensitivity_functions(
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

@testset "Test Gradients - OMIB; H; SourceBusVoltageChange; InterpolatingAdjoint" begin
    path = mktempdir()
    try
        omib_sys = build_system(PSIDTestSystems, "psid_test_omib")
        s_device = get_component(Source, omib_sys, "InfBus")
        s_change = SourceBusVoltageChange(1.0, s_device, :V_ref, 1.02)
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

        function plot_traces(δ, δ_gt)
            display(plot([scatter(; y = δ_gt), scatter(; y = δ)]))
        end
        EnzymeRules.inactive(::typeof(plot_traces), args...) = nothing
        function f_loss(states, δ_gt)
            #plot_traces(states[1], δ_gt)
            return sum(abs.(states[1] - δ_gt))
        end

        #GET SENSITIVITY FUNCTIONS 
        f_forward, f_grad, _ = get_sensitivity_functions(
            sim,
            [("generator-102-1", :Shaft, :H)],
            [("generator-102-1", :δ)],
            Rodas4(),
            f_loss;
            sensealg = InterpolatingAdjoint(; autojacvec = ReverseDiffVJP()),
            abstol = 1e-9,
            reltol = 1e-9,
            dtmax = 0.005,
            saveat = 0.005,
        )

        loss_zero = f_forward(p, δ_gt)
        loss_non_zero_1 = f_forward([3.2], δ_gt)
        loss_non_zero_2 = f_forward(p, δ_gt .* 2)
        @test loss_zero == 0.0
        @test loss_non_zero_1 == 0.36199910927656687
        @test loss_non_zero_2 == 172.66293171283323
        grad_zero = f_grad(p, δ_gt)
        grad_nonzero_1 = f_grad([3.14], δ_gt)
        grad_nonzero_2 = f_grad([3.15], δ_gt)
        @test isapprox(grad_zero[1], -1.0, atol = 1.0)
        @test isapprox(grad_nonzero_1[1], -8.0, atol = 1.0)
        @test isapprox(grad_nonzero_2[1], 8.0; atol = 1.0)
    finally
        @info("removing test files")
        rm(path; force = true, recursive = true)
    end
end

@testset "Test Gradients - OMIB; All; SourceBusVoltageChange" begin
    path = mktempdir()
    try
        omib_sys = build_system(PSIDTestSystems, "psid_test_omib")
        s_device = get_component(Source, omib_sys, "InfBus")
        s_change = SourceBusVoltageChange(1.0, s_device, :V_ref, 1.02)
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
        p = get_parameter_values(sim, :All)

        function plot_traces(δ, δ_gt)
            display(plot([scatter(; y = δ_gt), scatter(; y = δ)]))
        end
        EnzymeRules.inactive(::typeof(plot_traces), args...) = nothing
        function f_loss(states, δ_gt)
            #plot_traces(states[1], δ_gt)
            return sum(abs.(states[1] - δ_gt))
        end

        #GET SENSITIVITY FUNCTIONS 
        f_forward, f_grad, _ = get_sensitivity_functions(
            sim,
            :All,
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
        loss_non_zero_1 = f_forward(p .* 1.01, δ_gt)

        @test isapprox(loss_zero, 0.0, atol = 1e-9)
        @test isapprox(loss_non_zero_1, 1.49, atol = 1e-3)

        grad_zero = f_grad(p, δ_gt)
        @test isapprox(sum(grad_zero), 523.5340704294135, atol = 1.0)
    finally
        @info("removing test files")
        rm(path; force = true, recursive = true)
    end
end

@testset "Test Gradients - OMIB; H; PerturbState" begin
    path = mktempdir()
    try
        omib_sys = build_system(PSIDTestSystems, "psid_test_omib")
        s_device = get_component(Source, omib_sys, "InfBus")
        #p_branch = BranchImpedanceChange(1.0, Line,  "BUS 1-BUS 2-i_1", 2) BranchPerturbations not supported
        p_state = PerturbState(1.0, 5, 0.18)
        sim = Simulation!(
            MassMatrixModel,
            omib_sys,
            path,
            (0.0, 5.0),
            p_state,
        )
        #GET GROUND TRUTH DATA 
        execute!(sim, Rodas4(); abstol = 1e-9, reltol = 1e-9, dtmax = 0.005, saveat = 0.005)
        res = read_results(sim)
        t, δ_gt = get_state_series(res, ("generator-102-1", :δ); unique_timestamps = true)   #Avoid filtering of repeated timesteps 
        #GET PARAMETER VALUES 
        p = get_parameter_values(sim, [("generator-102-1", :Shaft, :H)])

        function plot_traces(δ, δ_gt)
            display(plot([scatter(; y = δ_gt), scatter(; y = δ)]))
        end
        EnzymeRules.inactive(::typeof(plot_traces), args...) = nothing
        function f_loss(states, δ_gt)
            #plot_traces(states[1], δ_gt)
            return sum(abs.(states[1] - δ_gt))
        end

        f_forward, f_grad, f_zygote_forward = get_sensitivity_functions(
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
        @test f_forward(p, δ_gt) == f_zygote_forward(p, δ_gt)
        @test f_forward([3.14], δ_gt) == f_zygote_forward([3.14], δ_gt)
        @test f_forward([3.15], δ_gt) == f_zygote_forward([3.15], δ_gt)

        @test f_grad(p, δ_gt) == Zygote.gradient(p -> f_zygote_forward(p, δ_gt), p)[1]

        _, _, f_zygote_forward = get_sensitivity_functions(
            sim,
            [("generator-102-1", :Shaft, :H)],
            [("generator-102-1", :δ)],
            Rodas4(),
            f_loss;
            sensealg = ReverseDiffAdjoint(),
            abstol = 1e-9,
            reltol = 1e-9,
            dtmax = 0.005,
            saveat = 0.005,
        )
        @test isapprox(
            Zygote.gradient(p -> f_zygote_forward(p, δ_gt), p)[1][1],
            78.6312919795057,
            atol = 1e-6,
        )
        @test isapprox(
            Zygote.gradient(p -> f_zygote_forward(p, δ_gt), [3.14])[1][1],
            53.94277234225185,
            atol = 1e-6,
        )
        @test isapprox(
            Zygote.gradient(p -> f_zygote_forward(p, δ_gt), [3.15])[1][1],
            78.48879776368774,
            atol = 1e-6,
        )
    finally
        @info("removing test files")
        rm(path; force = true, recursive = true)
    end
end

@testset "Test Optimization - OMIB; H" begin
    path = mktempdir()
    try
        omib_sys = build_system(PSIDTestSystems, "psid_test_omib")
        s_device = get_component(Source, omib_sys, "InfBus")
        s_change = SourceBusVoltageChange(1.0, s_device, :V_ref, 1.02)
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

        function plot_traces(δ, δ_gt)
            display(plot([scatter(; y = δ_gt), scatter(; y = δ)]))
        end
        EnzymeRules.inactive(::typeof(plot_traces), args...) = nothing
        function f_loss(states, δ_gt)
            #plot_traces(states[1], δ_gt)
            return sum(abs.(states[1] - δ_gt))
        end

        #GET SENSITIVITY FUNCTIONS 
        f_forward, f_grad, _ = get_sensitivity_functions(
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
        omib_sys = build_system(PSIDTestSystems, "psid_test_omib")
        s_device = get_component(Source, omib_sys, "InfBus")
        s_change = SourceBusVoltageChange(1.0, s_device, :V_ref, 1.02)
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
        function f_loss(states, δ_gt)
            #plot_traces(states[1], δ_gt)
            return sum(abs.(states[1] - δ_gt))
        end

        #GET SENSITIVITY FUNCTIONS 
        f_forward, f_grad, _ = get_sensitivity_functions(
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
        @test isapprox(grad_zero[1], 498.0, atol = 1.0)
        @test isapprox(grad_nonzero_1[1], 499.0, atol = 1.0)
        @test isapprox(grad_nonzero_2[1], -498.0; atol = 1.0)
    finally
        @info("removing test files")
        rm(path; force = true, recursive = true)
    end
end

function transform_power_load_to_constant_impedance(x::PowerLoad)
    l = StandardLoad(;
        name = get_name(x),
        available = get_available(x),
        bus = get_bus(x),
        base_power = get_base_power(x),
        constant_active_power = 0.0,
        constant_reactive_power = 0.0,
        impedance_active_power = get_active_power(x),
        impedance_reactive_power = get_reactive_power(x),
        current_active_power = 0.0,
        current_reactive_power = 0.0,
        max_constant_active_power = 0.0,
        max_constant_reactive_power = 0.0,
        max_impedance_active_power = get_max_active_power(x),
        max_impedance_reactive_power = get_max_reactive_power(x),
        max_current_active_power = 0.0,
        max_current_reactive_power = 0.0,
        services = get_services(x),
        dynamic_injector = get_dynamic_injector(x),
    )
    return l
end
using PlotlyJS
@time @testset "Test Gradients - 9 bus; Xd_p" begin
    path = mktempdir()
    try
        sys = build_system(PSIDTestSystems, "psid_test_ieee_9bus")
        for l in get_components(PowerLoad, sys)
            l_new = transform_power_load_to_constant_impedance(l)
            remove_component!(sys, l)
            add_component!(sys, l_new)
        end
        dyn_gen = get_component(DynamicGenerator, sys, "generator-3-1")
        get_machine(dyn_gen)
        sim = Simulation!(
            MassMatrixModel,
            sys,
            path,
            (0.0, 5.0),
            ControlReferenceChange(1.0, dyn_gen, :P_ref, 0.5),
        )

        #GET GROUND TRUTH DATA 
        execute!(sim, Rodas4(); abstol = 1e-6, reltol = 1e-6, dtmax = 0.05, saveat = 0.05)
        res = read_results(sim)
        t, δ_gt = get_state_series(res, ("generator-3-1", :δ))

        #GET PARAMETER VALUES 
        p = get_parameter_values(sim, [("generator-3-1", :Machine, :Xd_p)])

        function plot_traces(δ, δ_gt)
            display(plot([scatter(; y = δ_gt), scatter(; y = δ)]))
        end
        EnzymeRules.inactive(::typeof(plot_traces), args...) = nothing
        function f_loss(states, δ_gt)
            #plot_traces(states[1], δ_gt)
            return sum(abs.(states[1] - δ_gt))
        end

        #GET SENSITIVITY FUNCTIONS 
        f_forward, f_grad, _ = get_sensitivity_functions(
            sim,
            [("generator-3-1", :Machine, :Xd_p)],
            [("generator-3-1", :δ)],
            Rodas4(),
            f_loss;
            sensealg = InterpolatingAdjoint(; autojacvec = ReverseDiffVJP()),
            abstol = 1e-6,
            reltol = 1e-6,
            dtmax = 0.05,
            saveat = 0.05,
        )

        loss_zero = f_forward(p, δ_gt)
        loss_non_zero_1 = f_forward(p * 1.01, δ_gt)
        loss_non_zero_2 = f_forward(p * 0.99, δ_gt)
        @test isapprox(loss_zero, 0.0, atol = 2e-9)
        @test loss_non_zero_1 != 0.0
        @test loss_non_zero_2 != 0.0
        grad_zero = f_grad(p, δ_gt)
        grad_nonzero_1 = f_grad(p * 1.01, δ_gt)
        grad_nonzero_2 = f_grad(p * 0.99, δ_gt)
        @test isapprox(grad_zero[1], 0.8876855633591412, atol = 1e-6)
        @test isapprox(grad_nonzero_1[1], 0.5946989856464833, atol = 1.0)
        @test isapprox(grad_nonzero_2[1], -0.9432527319604077; atol = 1.0)
    finally
        @info("removing test files")
        rm(path; force = true, recursive = true)
    end
end

#= 
#BELOW HERE: TESTS FOR SENSITIVITY OF DDES

function add_degov_to_omib!(omib_sys)
    gen = get_component(ThermalStandard, omib_sys, "generator-102-1")
    dyn_gen = get_component(DynamicGenerator, omib_sys, "generator-102-1")
    new_gov = PSY.DEGOV(;
        T1 = 0.0,
        T2 = 0.0,
        T3 = 0.0,
        K = 18.0,
        T4 = 12.0,
        T5 = 5.0,
        T6 = 0.2,
        Td = 0.5,
        P_ref = 0.0,
    )
    dyn_gen_new = DynamicGenerator(;
        name = get_name(dyn_gen),
        ω_ref = get_ω_ref(dyn_gen),
        machine = get_machine(dyn_gen),
        shaft = get_shaft(dyn_gen),
        avr = get_avr(dyn_gen),
        prime_mover = new_gov,
        pss = get_pss(dyn_gen),
        base_power = get_base_power(dyn_gen),
    )
    remove_component!(omib_sys, dyn_gen)
    add_component!(omib_sys, dyn_gen_new, gen)
end

 @testset "Test Gradients - OMIB; H; Delays" begin
    path = mktempdir()
    try
        #EnzymeRules.inactive(::typeof(Base.hasproperty), args...) = nothing # To allow wrap_sol to work?  --> causes gradient to b zero; bad idea to go marking functions invalid without clear reason. 
        omib_sys = build_system(PSIDTestSystems, "psid_test_omib")
        s_device = get_component(Source, omib_sys, "InfBus")
        s_change = SourceBusVoltageChange(1.0, s_device, :V_ref, 1.02)
        add_degov_to_omib!(omib_sys)
        sim = Simulation!(
            MassMatrixModel,
            omib_sys,
            pwd(),
            (0.0, 5.0), #Inside of the delay time... 
            s_change,
        )

        #GET GROUND TRUTH DATA 
        execute!(
            sim,
            MethodOfSteps(Rodas5(; autodiff = false));
            abstol = 1e-9,
            reltol = 1e-9,
            dtmax = 0.005,
            saveat = 0.005,
        )
        res = read_results(sim)
        t, δ_gt = get_state_series(res, ("generator-102-1", :δ))

        #GET PARAMETER VALUES 
        p = get_parameter_values(sim, [("generator-102-1", :Shaft, :H)])

        function plot_traces(δ, δ_gt)
            display(plot([scatter(; y = δ_gt), scatter(; y = δ)]))
        end
        EnzymeRules.inactive(::typeof(plot_traces), args...) = nothing
        function f_loss(states, δ_gt)
            #plot_traces(δ, δ_gt)
            return sum(abs.(states[1] - δ_gt))
        end

        #GET SENSITIVITY FUNCTIONS 
        f_forward, f_grad, _ = get_sensitivity_functions(
            sim,
            [("generator-102-1", :Shaft, :H)],
            [("generator-102-1", :δ)],
            MethodOfSteps(Rodas4()), #; autodiff = false
            #MethodOfSteps(QNDF(;autodiff=true)),
            f_loss;
            sensealg = ForwardDiffSensitivity(),
            abstol = 1e-6, # 1e-6
            reltol = 1e-6, # 1e-6 
            # dtmax = 0.005,
             saveat = 0.005,
        )
       # @error length(δ_gt)
        loss_zero = f_forward(p, δ_gt)
        loss_non_zero_1 = f_forward([5.2], δ_gt)
        loss_non_zero_2 = f_forward(p, δ_gt .* 2)
        @test isapprox(loss_zero, 0.0, atol = 3e-3)
        @test loss_non_zero_1 != 0.0
        @test loss_non_zero_2 != 0.0
        grad_zero = f_grad(p, δ_gt)
        @error grad_zero
        @test isapprox(grad_zero[1], 0.0, atol = 1e-12)
    finally
        @info("removing test files")
        rm(path; force = true, recursive = true)
    end
end

#Potential Rank Deficient Matrix Detected: Erorrs... 
@testset "Test Gradients - OMIB; Xd_p; delays" begin
    path = mktempdir()
    try
        omib_sys = build_system(PSIDTestSystems, "psid_test_omib")
        add_degov_to_omib!(omib_sys)
        s_device = get_component(Source, omib_sys, "InfBus")
        s_change = SourceBusVoltageChange(1.0, s_device, :V_ref, 1.02)
        sim = Simulation!(
            MassMatrixModel,
            omib_sys,
            path,
            (0.0, 0.05),
            #s_change,
        )

        #GET GROUND TRUTH DATA 
        execute!(
            sim,
            MethodOfSteps(Rodas4());
            abstol = 1e-9,
            reltol = 1e-9,
            #dtmax = 0.005,
            saveat = 0.005,
        )
        res = read_results(sim)
        t, δ_gt = get_state_series(res, ("generator-102-1", :δ))

        #GET PARAMETER VALUES 
        p = get_parameter_values(sim, [("generator-102-1", :Machine, :Xd_p)])

        function plot_traces(δ, δ_gt)
            display(plot([scatter(; y = δ_gt), scatter(; y = δ)]))
        end
        EnzymeRules.inactive(::typeof(plot_traces), args...) = nothing
        function f_loss(states, δ_gt)
            #plot_traces(δ, δ_gt)
            return sum(abs.(states[1] - δ_gt))
        end

        #GET SENSITIVITY FUNCTIONS 
        f_forward, f_grad = get_sensitivity_functions(
            sim,
            [("generator-102-1", :Machine, :Xd_p)],
            [("generator-102-1", :δ)],
            MethodOfSteps(Rodas4()),
            f_loss;
            sensealg = ForwardDiffSensitivity(),
            abstol = 1e-6,
            reltol = 1e-6,
           # dtmax = 0.005,
            saveat = 0.005,
        )

        loss_zero = f_forward(p, δ_gt)
        #loss_non_zero_1 = f_forward(p * 1.01, δ_gt)
        #loss_non_zero_2 = f_forward(p * 0.99, δ_gt)
        #@test isapprox(loss_zero, 0.0, atol = 1e-9)
        #@test loss_non_zero_1 != 0.0
        #@test loss_non_zero_2 != 0.0
        #grad_zero = f_grad(p, δ_gt)
        #grad_nonzero_1 = f_grad(p * 1.01, δ_gt)
        #grad_nonzero_2 = f_grad(p * 0.99, δ_gt)
        #@test isapprox(grad_zero[1], 497, atol = 1.0)
        #@test isapprox(grad_nonzero_1[1], 499.0, atol = 1.0)
        #@test isapprox(grad_nonzero_2[1], -498.0; atol = 1.0)
    finally
        @info("removing test files")
        rm(path; force = true, recursive = true)
    end
end

 =#
