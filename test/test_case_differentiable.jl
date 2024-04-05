"""
TODO: Make a test for the simplest case of taking gradient wrt parameters (add SciMLSensitivity and Zygote to test dependencies)
-Does this test belong in PSID?
-OMIB system.
-Solve with H=3, save ground truth data.
-Build function that sets H=2, builds simulation, executes, and calculates loss relative to prior result.
-Take gradient with Zygote. 
"""

#Write a test which differentiates through a call to _simple_execute for the OMIB system.
#IF we can actually calculate a gradient, compare to a gradient calculated offline from "scratch" system. 
#If we get this working, I'd say we have the "core" changes done. 
#Create an issue (show the use case in pseudo-code) and a PR in github. 

##################################################
############### LOAD DATA ########################
##################################################

using Revise
using PowerSimulationsDynamics
using PowerSystems
using Test
using NLsolve
using SciMLBase
using Sundials
using OrdinaryDiffEq
using DelimitedFiles
using DataFrames
using InfrastructureSystems
using PowerSystemCaseBuilder
using PowerFlows
using PowerNetworkMatrices
import LinearAlgebra
using Logging
using Zygote
using SciMLSensitivity
using PlotlyJS
import ChainRulesCore
import TimerOutputs
#import Aqua
#Aqua.test_unbound_args(PowerSimulationsDynamics)
#Aqua.test_undefined_exports(PowerSimulationsDynamics)
#Aqua.test_ambiguities(PowerSimulationsDynamics)
#Aqua.test_stale_deps(PowerSimulationsDynamics)
#Aqua.test_deps_compat(PowerSimulationsDynamics)
const IS = InfrastructureSystems
const PSY = PowerSystems
const PSID = PowerSimulationsDynamics
const PNM = PowerNetworkMatrices
const PSB = PowerSystemCaseBuilder
const CRC = ChainRulesCore
Zygote.refresh()    #similar to restart of REPL after changes? Seemingly doesn't work for changes to PSID source sometimes (try first, but if unexpected also restart REPL)
tspan = (0.0, 5.0)
t_perturb = 1.0
Pref_perturb = 1.0
println()
# omib system, add dynamic lines
sys = build_system(PSIDTestSystems, "psid_test_omib")
for b in get_components(Line, sys)
    dyn_branch = DynamicBranch(b)
    add_component!(sys, dyn_branch)
end

#SOLVE SIMULATION USING PSID FOR TWO VALUES OF H
dyngen = get_component(DynamicGenerator, sys, "generator-102-1")
PSY.set_H!(PSY.get_shaft(dyngen), 3.0)
pert = ControlReferenceChange(t_perturb, dyngen, :P_ref, Pref_perturb)
sim = Simulation(MassMatrixModel, sys, pwd(), tspan, pert)
execute!(sim, Rodas5(); saveat = 0.01)
res = read_results(sim)
t_gt, δ_gt = get_state_series(res, ("generator-102-1", :δ); dt = 0.01)

@test typeof(t_gt) == Vector{Float64}
@test typeof(δ_gt) == Vector{Float64}

#TODO - write a test which loops through and tries sensitivity ooptions. 
# Zygote less restrictive on types than ReverseDiff? 
# Prefer ZygoteVJP() and EnzymeVJP()



#ReverseDiffVJP and EnzymeVJP only compatible with Hybrid DEs (DEs with callbacks)

#However, BacksolveAdjoint is prone to instabilities whenever the Lipschitz constant is sufficiently large, 
#like in stiff equations, PDE discretizations, and many other contexts, so it is not used by default.

sensealgs = (
    InterpolatingAdjoint(checkpointing=true, autojacvec = ReverseDiffVJP(true)),
    #InterpolatingAdjoint(autojacvec = EnzymeVJP()),
    #InterpolatingAdjoint(autojacvec = ZygoteVJP()),
    #BacksolveAdjoint(autojacvec = ReverseDiffVJP(true)),
    #BacksolveAdjoint(autojacvec = ZygoteVJP()),
    #BacksolveAdjoint(autojacvec = ReverseDiffVJP(false)),
    #BacksolveAdjoint(autojacvec = TrackerVJP()),
    #QuadratureAdjoint(autojacvec = ReverseDiffVJP(true)),
    #TrackerAdjoint()
)
 



for sensealg in sensealgs
    dyngen = collect(get_components(DynamicGenerator, sys))[1]
    function get_v1(_p)
        #buf = Zygote.Buffer(sim.problem.p)  
        #for i in eachindex(buf)
        #    if i == 4
        #        buf[i] = θ[1]
        #    else 
        #        buf[i] = sim.problem.p[i]
        #    end 
        #end 
        #new_p = copy(buf)

        set_H!(get_shaft(dyngen), 2.0)
        PSID._build_inputs!(sim)
        prob_ = remake(sim.problem; p = _p)
        sim.problem = prob_

        #saveat problematic? 
        PSID._simple_execute!(sim, Rodas5(autodiff=false); saveat=0.1, sensealg = sensealg) 
        _, v101 = PSID.get_voltage_magnitude_series(sim.results, 101)
        return v101
    end 
    function loss(_p)
        v101_pred = get_v1(_p)
        return sum(abs.(v101_pred - data_gt))
    end 
    
    sim = Simulation(MassMatrixModel, sys, pwd(), tspan, pert)  #Need this to reset P_ref 
    p_init = sim.problem.p
    data_gt = get_v1(p_init)



        #display(plot([
        #    scatter(y=aa)
        #])) 
   
        #=     ts = collect(range(0; stop = sim.results.solution.t[end], step = 0.01))
            state = sim.results.solution(collect(ts); idxs = 7)
            δ_pred = state.u
             display(plot([
                scatter(x = ts, y=δ_pred)
                scatter(x = t_gt, y=δ_gt)
            ])) 
            loss = sum(abs.(δ_pred .- δ_gt)) 
            return loss =#


    println(length(p_init))
    loss(p_init)
    println("loss with original parameters   ", loss(p_init .* 1.1 ))
    #println("loss with parameters + 1%   ", loss(p_init*1.01))
    #println("loss with parameters - 1%   ", loss(p_init*0.99))
    #println(loss([2.0]))
    #println(loss([2.0]))
    #println(loss([3.0]))
    #println(loss([3.0]))
    #println(loss([4.0]))
    #println(loss([4.0]))
    ###@test isapprox(loss([3.0]), 0.0)

    #sim = Simulation(MassMatrixModel, sys, pwd(), tspan, pert)  #Need this to reset P_ref 
    #@test isapprox(loss([2.0]), 39.992118700338644)

    #sim = Simulation(MassMatrixModel, sys, pwd(), tspan, pert)  #Need this to reset P_ref 
    #@test isapprox(loss([4.0]), 40.742272618402)

    #CRC.@ignore_derivatives @error SciMLSensitivity.automatic_sensealg_choice(sim.problem, sim.problem.u0, sim.problem.p, true)      #Deterimine default sensealg
  
    try 
        xx = @timed Zygote.gradient((p) -> loss(p), p_init)
        println("Success: $sensealg   gradient $(xx.value) took $(xx.time)s")
    catch e
         @error "Failure: $sensealg"
        Base.showerror(stdout, e)
        Base.show_backtrace(stdout, Base.catch_backtrace())
        println()
        println()
    end 
end 

#Gradient wrt H at H = 4 ([-5.931909504995866],) - these should be approx
#Gradient wrt H at H = 2 ([23.147017043034136 


