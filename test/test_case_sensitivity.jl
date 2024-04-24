
"""
Case 1:
This case study defines a classical machine against an infinite bus. Sensitivitiy
Analysis is performed for the inertia of the generator. The fault is a change in the 
source voltage.
"""
##################################################
############### LOAD DATA ######################                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               ##
##################################################

omib_sys = build_system(PSIDTestSystems, "psid_test_omib")

##################################################
############### SOLVE PROBLEM ####################
##################################################

s_device = get_component(Source, omib_sys, "InfBus")
s_change = SourceBusVoltageChange(1.0, s_device, :V_ref, 1.02)
#using PlotlyJS  #for debug only 


#NOTES ON SENSITIVITY ALGORITHMS FROM SCIMLSENSITIVITY               
#ReverseDiffVJP and EnzymeVJP only options compatible with Hybrid DEs (DEs with callbacks)
#BacksolveAdjoint prone to instabilities whenever the Lipschitz constant is sufficiently large (stiff equations, PDE discretizations, and many other contexts) 
#For ForwardDiffSensitivity, convert_tspan=true is needed for hybrid equations. 
@testset "Test Gradients - Mass Matrix no delays" begin
    path = mktempdir()
    try
        sim = Simulation!(
            MassMatrixModel,
            omib_sys,
            path,
            (0.0, 5.0),
            s_change,
        )
        H_gt =
            get_H(get_shaft(get_component(DynamicGenerator, omib_sys, "generator-102-1")))
        execute!(sim, Rodas5(); dtmax = 0.005, saveat = 0.005)
        res = read_results(sim)
        t, δ_gt = get_state_series(res, ("generator-102-1", :δ))
        for solver in [FBDF(), Rodas5(), QNDF()]
            for tol in [1e-6, 1e-9, 1e-12]
                function f(sim)
                    execute!(
                        sim,
                        solver;
                        abstol = tol,
                        reltol = tol,
                        dtmax = 0.005,
                        saveat = 0.005,
                    )
                    res = read_results(sim)
                    t, δ = get_state_series(res, ("generator-102-1", :δ))
                    #display(plot(scatter(x=t, y = δ)))
                    return sum(abs.(δ - δ_gt))
                end
                g = PSID.get_parameter_sensitivity_function!(
                    sim,
                    [("generator-102-1", SingleMass, :H)],
                    f,
                )
                #p = PSID.get_parameter_sensitivity_values(sim, [("generator-102-1", SingleMass, :H)])
                @test isapprox(Zygote.gradient(g, [3.14])[1][1], -8.0, atol = 1.0)
                @test isapprox(Zygote.gradient(g, [3.15])[1][1], 8.0, atol = 1.0)
            end
        end
    finally
        @info("removing test files")
        rm(path; force = true, recursive = true)
    end
end

@testset "Test Parameter Optimization Problem - no delays" begin
    path = mktempdir()
    try
        sim = Simulation!(
            MassMatrixModel,
            omib_sys,
            path,
            (0.0, 5.0),
            s_change,
        )
        H_gt =
            get_H(get_shaft(get_component(DynamicGenerator, omib_sys, "generator-102-1")))
        execute!(sim, Rodas5(; autodiff = true); dtmax = 0.005, saveat = 0.005)
        res = read_results(sim)
        t, δ_gt = get_state_series(res, ("generator-102-1", :δ))

        function f(sim)
            execute!(sim, FBDF(; autodiff = true); dtmax = 0.005, saveat = 0.005)
            res = read_results(sim)
            t, δ = get_state_series(res, ("generator-102-1", :δ))
            #display(plot(scatter(x=t, y = δ)))
            return sum(abs.(δ - δ_gt))
        end
        g = PSID.get_parameter_sensitivity_function!(
            sim,
            [("generator-102-1", SingleMass, :H)],
            f,
        )
        p = PSID.get_parameter_sensitivity_values(
            sim,
            [("generator-102-1", SingleMass, :H)],
        )
        #H_values = [] 
        #loss_values = []
        function callback(u, l)
            #push!(H_values, u.u[1])
            #push!(loss_values, l)
            return false
        end
        optfun = OptimizationFunction((u, _) -> g(u), Optimization.AutoZygote())
        optprob = OptimizationProblem(optfun, [3.14])
        sol = Optimization.solve(
            optprob,
            OptimizationOptimisers.Adam(0.002);
            callback = callback,
            maxiters = 3,
        )
        @test sol.u[1] == 3.144000187737475
        #display(plot(scatter(y=H_values)))
        #display(plot(scatter(y=loss_values)))
    finally
        @info("removing test files")
        rm(path; force = true, recursive = true)
    end
end

@testset "Test Gradients - Mass Matrix with delays" begin
    path = mktempdir()
    try
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

        sim = Simulation!(
            MassMatrixModel,
            omib_sys,
            path,
            (0.0, 5.0),
            s_change,
        )
        H_gt =
            get_H(get_shaft(get_component(DynamicGenerator, omib_sys, "generator-102-1")))
        execute!(
            sim,
            MethodOfSteps(Rodas5(; autodiff = false));
            dtmax = 0.005,
            saveat = 0.005,
        )
        res = read_results(sim)
        t, δ_gt = get_state_series(res, ("generator-102-1", :δ))
        for solver in [
            MethodOfSteps(Rodas5(; autodiff = true)),
            MethodOfSteps(QNDF(; autodiff = true)),
        ]
            for tol in [1e-6]
                function f(sim)
                    execute!(
                        sim,
                        solver;
                        abstol = tol,
                        reltol = tol,
                        dtmax = 0.005,
                        saveat = 0.005,
                    )
                    res = read_results(sim)
                    t, δ = get_state_series(res, ("generator-102-1", :δ))
                    #display(plot(scatter(x=t, y = δ)))
                    return sum(abs.(δ - δ_gt))
                end
                g = PSID.get_parameter_sensitivity_function!(
                    sim,
                    [("generator-102-1", SingleMass, :H)],
                    f,
                )
                #p = PSID.get_parameter_sensitivity_values(sim, [("generator-102-1", SingleMass, :H)])
                #display(Zygote.gradient(g, [3.14]))
                @test isapprox(Zygote.gradient(g, [3.14])[1][1], -10.0, atol = 1.0)
                @test isapprox(Zygote.gradient(g, [3.15])[1][1], 10.0, atol = 1.0)
            end
        end
    finally
        @info("removing test files")
        rm(path; force = true, recursive = true)
    end
end

@testset "Test Parameter Optimization Problem - with delays" begin
    path = mktempdir()
    try
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
        sim = Simulation!(
            MassMatrixModel,
            omib_sys,
            path,
            (0.0, 5.0),
            s_change,
        )
        H_gt =
            get_H(get_shaft(get_component(DynamicGenerator, omib_sys, "generator-102-1")))
        execute!(
            sim,
            MethodOfSteps(Rodas5(; autodiff = true));
            dtmax = 0.005,
            saveat = 0.005,
        )
        res = read_results(sim)
        t, δ_gt = get_state_series(res, ("generator-102-1", :δ))

        function f(sim)
            execute!(
                sim,
                MethodOfSteps(Rodas5(; autodiff = true));
                abstol = 1e-6,
                reltol = 1e-6,
                dtmax = 0.005,
                saveat = 0.005,
            )
            res = read_results(sim)
            t, δ = get_state_series(res, ("generator-102-1", :δ))
            #display(plot(scatter(x=t, y = δ)))
            return sum(abs.(δ - δ_gt))
        end
        g = PSID.get_parameter_sensitivity_function!(
            sim,
            [("generator-102-1", SingleMass, :H)],
            f,
        )
        p = PSID.get_parameter_sensitivity_values(
            sim,
            [("generator-102-1", SingleMass, :H)],
        )
        #H_values = [] 
        #loss_values = []
        function callback(u, l)
            #push!(H_values, u.u[1])
            #push!(loss_values, l)
            return false
        end
        optfun = OptimizationFunction((u, _) -> g(u), Optimization.AutoZygote())
        optprob = OptimizationProblem(optfun, [3.14])
        sol = Optimization.solve(
            optprob,
            OptimizationOptimisers.Adam(0.002);
            callback = callback,
            maxiters = 3,
        )
        @test sol.u[1] == 3.1440015515797763
        #display(plot(scatter(y=H_values)))
        #display(plot(scatter(y=loss_values)))
    finally
        @info("removing test files")
        rm(path; force = true, recursive = true)
    end
end
