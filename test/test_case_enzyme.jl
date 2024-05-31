
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

sim = Simulation(
    MassMatrixModel,
    omib_sys,
    pwd(),
    (0.0, 5.0),
    s_change,
)

#JULIA CRASHES  - the problem has the Jacobian... 
#= function f(p_in, sim )
    prob = sim.problem
    prob_new = remake(prob, u0 = vcat(sim.problem.u0[1:end-1], p_in))
    #display(prob_new)
    sol = solve(prob_new, Rodas4(), abstol =1e-9, reltol=1e-9)
    u1 = [u[1] for u in sol.u]
    #u6 = [u[6] for u in sol.u]
    return sum(abs.(u1))
end 
p = [1.0000001]
f(p, sim)
dp = make_zero(p)
dsim = make_zero(sim)
#Would be huge if this works! 
=#
#= prob = sim.problem
function f(p_in, prob_in )
    prob_new = remake(prob_in, u0 = vcat(prob_in.u0[1:end-1], p_in))
    sum(solve(prob_new, Rodas4(), abstol =1e-9, reltol=1e-9))
end
p = [1.0000001]
f(p, prob)
dp_zygote = Zygote.gradient((p)->f(p, prob), p) =#

#THIS EXACT THING WORKS IF THE PROBLEM DOESN"T COME FROM A SIMULATION... 
prob = sim.problem
ode_f = prob.f
new_f = ODEFunction(
    ode_f.f;
    mass_matrix = ode_f.mass_matrix,
    jac = nothing,
    tgrad = ode_f.tgrad,
)
new_prob = ODEProblem(new_f, prob.u0, prob.tspan, prob.p)
function f(p_in, prob_in)
    prob_new = remake(prob_in; u0 = vcat(prob_in.u0[1:(end - 1)], p_in))
    sum(solve(prob_new, Rodas4(); abstol = 1e-9, reltol = 1e-9))
end
p = [1.0000002]
f(p, new_prob)
dp = make_zero(p)
dprob = make_zero(new_prob)
Enzyme.autodiff(Reverse, f, Active, Duplicated(p, dp), Duplicated(new_prob, dprob))  #Would be huge if this works! 
dp_zygote2 = Zygote.gradient((p) -> f(p, new_prob), p)[1]
dp
dp_zygote2

##Now pass in the full sim...
sim = Simulation(
    MassMatrixModel,
    omib_sys,
    pwd(),
    (0.0, 5.0),
    s_change,
)
sim.problem.f.jac
prob = sim.problem
ode_f = prob.f
new_f = ODEFunction(
    ode_f.f;
    mass_matrix = ode_f.mass_matrix,
    jac = nothing,
    tgrad = ode_f.tgrad,
)
new_prob = ODEProblem(new_f, prob.u0, prob.tspan, prob.p)
sim.problem = new_prob
sim.problem.f.jac

#= #JULIA CRASHES... 
#f2 = sim.problem.f.f #system model... Just a function. 
function f(p_in, sim_in)
    prob = sim_in.problem
    #ode_f = prob.f
    #new_f  = ODEFunction{true}(f2)#; mass_matrix = ode_f.mass_matrix, jac = nothing, tgrad = ode_f.tgrad)
    #new_prob = ODEProblem{true}(new_f, prob.u0, prob.tspan, prob.p)
    prob_new = remake(prob, u0 = vcat(prob.u0[1:end-1], p_in))
    sum(solve(prob_new, Rodas4(), abstol =1e-9, reltol=1e-9))
end
p = [1.0000002]
f(p, sim)
dp = make_zero(p)
dsim= make_zero(sim)
Enzyme.autodiff(Reverse, f, Active, Duplicated(p, dp), Duplicated(sim, dsim))  #Would be huge if this works! 
dp_zygote2 = Zygote.gradient((p)->f(p, sim), p)[1]
dp
dp_zygote2

## =#
#PASS ONLY PROBLEM AND SIM INPUTS... In case it is some other part of Sim that causes the problem . 
#=
 fds 
=#

#=  
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
        execute!(sim, Rodas5(); abstol = 1e-9, reltol = 1e-9, dtmax = 0.005, saveat = 0.005)
        res = read_results(sim)
        t, δ_gt = get_state_series(res, ("generator-102-1", :δ))
        for solver in [FBDF(), Rodas5(), QNDF()]
            for tol in [1e-6, 1e-9]
                function f(sim, δ_gt)
                    execute!(
                        sim,
                        solver;
                        sensealg = ForwardDiffSensitivity(),
                        abstol = tol,
                        reltol = tol,
                        dtmax = 0.005,
                        saveat = 0.005,
                    )
                    res = read_results(sim)
                    t, δ = get_state_series(res, ("generator-102-1", :δ))
                    display(plot(scatter(; x = t, y = δ)))
                    return sum(abs.(δ - δ_gt))
                end
                p = get_parameter_values(sim, [("generator-102-1", :Shaft, :H)])
                f_forward =
                    get_forward_function(f, sim, δ_gt, [("generator-102-1", :Shaft, :H)])
                @error p
                @error f_forward([3.148])
                @error f_forward([3.0])
                #@error Enzyme.make_zero(sim)
                f_grad =
                    get_gradient_function(f, sim, δ_gt, [("generator-102-1", :Shaft, :H)])
                f_grad([3.148])
                #f_grad = get_tradient_function
                #p = PSID.get_parameter_sensitivity_values(sim, [("generator-102-1", SingleMass, :H)])
                #@error Zygote.gradient(g, [3.15])[1][1]
                #=                 @test isapprox(
                                    Zygote.gradient((p) -> g(p, δ_gt), [3.14])[1][1],
                                    -8.0,
                                    atol = 1.0,
                                )
                                @test isapprox(
                                    Zygote.gradient((p) -> g(p, δ_gt), [3.15])[1][1],
                                    8.0,
                                    atol = 1.0,
                                ) =#
            end
        end
    finally
        @info("removing test files")
        rm(path; force = true, recursive = true)
    end
end
 =#
