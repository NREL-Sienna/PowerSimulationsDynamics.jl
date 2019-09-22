##Construct Dynamical System
OMIB4 = DynamicSystem(nodes_OMIB, #2bus Network
                      branch_OMIB, #Network data
                      [Gen2AVR], #4StateGen with AVR
                      [inf_gen], #Inf bus and load
                      100.0, #MVABase
                      60.0) #f_sys

##Initial conditions
dx0 = zeros(OMIB4.counts[:total_rows])
x0 = [1.1, #v_1R
    1.0, #v_2R
    0.0, #v_1I
    0.0, #v_2I
    1.05, #eq_p
    0.10, #ed_pp
    0.0, #δ
    1.0, #ω,
    1.05] #V_f
u0 = [0.4] #Pd
diff_vars = OMIB4.DAE_vector;
tspan = (0.0, 15.0);

##Find equilibrium point
inif! = (out,x) -> system_model!(out, dx0 ,x, (u0,OMIB4), 0.0)
sys_solve4 = nlsolve(inif!, x0)
x0_init = sys_solve4.zero

##Construct DAE problem
prob = DiffEqBase.DAEProblem(system_model!, dx0, x0_init, tspan,
                            (u0, OMIB4), differential_vars = diff_vars)

##Solve and plot DAE
sol = solve(prob,IDA(init_all = :false))
fig1 = plot(sol)

##Construct Callback
tstop = [1.0]
cb = DiffEqBase.DiscreteCallback(change_t_one, step_change3!)

##Solve problem with Callback (increase in Pd at t=1.0 from 0.4 to 0.5)
sol2 = solve(prob, IDA(init_all = :false), callback=cb, tstops=tstop)
fig2 = plot(sol2)

##Save Figure
savefig(fig2, "plots/sims/OMIB_4states_fixedAVR.png")

##Return turbine governor power to normal state
OMIB4.dyn_injections[1].P_ref = 0.4
