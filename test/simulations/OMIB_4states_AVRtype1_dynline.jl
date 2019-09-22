##Construct Dynamical System
OMIB5 = DynamicSystem(nodes_OMIB, #2bus Network
                      Vector{PSY.Line}(), #Network data
                      [Gen3AVR], #4StateGen with AVR
                      [inf_gen], #Inf bus and load
                      100.0, #MVABase
                      60.0,
                      Dynbranch_OMIB) #f_sys

##Initial conditions
dx0 = zeros(14)
x0 = [  1.049992525397527
1.0027262415732132
1.1049625442966575e-5
0.1602810688176053
0.9535007165758107
0.26150065989229665
0.4365934439055868
1.0
0.9076610164638759
0.908893689474798
-0.0010891932197566508
1.0154555315526261
-2.2099247968832048
-1.5226402972437107  ]
u0 = [0.4] #Pd
diff_vars = OMIB5.DAE_vector;
tspan = (0.0, 30.0)

##Find equilibrium point
inif! = (out, x) -> system_model!(out, dx0 , x, (u0, OMIB5), 0.0)
sys_solve5 = nlsolve(inif!, x0)
x0_init_dyn = sys_solve5.zero

##Construct DAE problem
prob = DiffEqBase.DAEProblem(system_model!, dx0, x0_init_dyn, tspan,
                            (u0, OMIB5), differential_vars = diff_vars)

##Solve and plot DAE
sol = solve(prob, IDA())
fig1 = plot(sol)

##Construct Callback
tstop = [1.0]
cb = DiffEqBase.DiscreteCallback(change_t_one, step_change3!)

##Solve problem with Callback (increase in Pd at t=1.0 from 0.4 to 0.5)
sol2 = solve(prob, IDA(init_all = :false), callback=cb, tstops=tstop)
fig2 = plot(sol2)

##Save Figure
savefig(fig2, "plots/sims/OMIB_4states_AVRtype1_dyn.png")

##Return turbine governor power to normal state
OMIB5.dyn_injections[1].P_ref = 0.4
