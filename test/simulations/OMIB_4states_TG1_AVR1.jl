##Construct Dynamical System
OMIB6 = DynamicSystem(nodes_OMIB, #2bus Network
                      branch_OMIB, #Network data
                      [Gen4AVRTG], #4StateGen with AVR
                      [inf_gen], #Inf bus and load
                      100.0, #MVABase
                      60.0) #f_sys

##Initial conditions
dx0 = zeros(OMIB6.counts[:total_rows])
x0 = [1.1, #v_1R
    1.0, #v_2R
    0.0, #v_1I
    0.0, #v_2I
    1.05, #eq_p
    0.10, #ed_pp
    0.0, #δ
    1.0, #ω,
    1.05, #Vf AVR
    0.5, #Vr1 AVR
    0.0, #Vr2 AVR
    1.05, #Vm AVR
    0.4, #x_g1 TG
    0.4, #x_g2 TG
    0.4] #x_g3 TG
u0 = [0.4] #Pd
diff_vars = OMIB6.DAE_vector;
tspan = (0.0, 200.0)

##Find equilibrium point
inif! = (out,x) -> system_model!(out, dx0 ,x, (u0,OMIB6), 0.0)
sys_solve6 = nlsolve(inif!, x0)
x0_init = sys_solve6.zero

##Construct DAE problem
prob = DiffEqBase.DAEProblem(system_model!, dx0, x0_init, tspan,
                          (u0, OMIB6), differential_vars = diff_vars)

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
savefig(fig2, "plots/sims/OMIB_4states_TG1_AVR1.png")

##Return turbine governor power to normal state
OMIB6.dyn_injections[1].P_ref = 0.4
