##Construct Dynamical System
OMIB5 = DynamicSystem(nodes_OMIB, #2bus Network
                      branch_OMIB, #Network data
                      [Gen3AVR], #4StateGen with AVR
                      [inf_gen], #Inf bus and load
                      100.0, #MVABase
                      60.0) #f_sys

##Initial conditions
dx0 = zeros(get_total_rows(OMIB5))
x0 = [1.1, #v_1R
    1.0, #v_2R
    0.0, #v_1I
    0.0, #v_2I
    1.05, #eq_p
    0.10, #ed_pp
    0.0, #δ
    1.0, #ω,
    1.05, #Vf
    0.5, #Vr1
    0.0, #Vr2
    1.05] #Vm
u0 = [0.4] #Pd
diff_vars = OMIB5.DAE_vector;
tspan = (0.0, 30.0)

##Find equilibrium point
inif! = (out,x) -> system_model!(out, dx0 ,x, (u0,OMIB5), 0.0)
sys_solve5 = nlsolve(inif!, x0)
x0_init = sys_solve5.zero

##Construct DAE problem
prob = DiffEqBase.DAEProblem(system_model!, dx0, x0_init, tspan,
                            (u0, OMIB5), differential_vars = diff_vars)

##Solve and plot DAE
sol = solve(prob,IDA())
fig1 = plot(sol)

##Construct Callback
tstop = [1.0]
cb = DiffEqBase.DiscreteCallback(change_t_one, step_change3!)

##Solve problem with Callback (increase in Pd at t=1.0 from 0.4 to 0.5)
sol2 = solve(prob, IDA(init_all = :false), callback=cb, tstops=tstop)
fig2 = plot(sol2)

##Save Figure
savefig(fig2, "plots/sims/OMIB_4states_AVRtype1.png")

##Return turbine governor power to normal state
OMIB5.dyn_injections[1].P_ref = 0.4
