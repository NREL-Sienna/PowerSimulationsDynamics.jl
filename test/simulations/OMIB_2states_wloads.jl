##Construct Dynamical System
OMIB2 = DynamicSystem(nodes_OMIB, #2bus Network
                      branch_OMIB, #Network data
                      [Gen1AVRnoAVR], #2StateGen without AVR
                      [inf_gen, loads[1]], #Inf bus and load
                      100.0, #MVABase
                      60.0) #f_sys

##Initial conditions
dx0 = zeros(OMIB2.counts[:total_rows])
x0 = [1.1, 1.0, 0.0, 0.0, 0.0, 1.0]
u0 = [0.4] #Pd
diff_vars = OMIB2.DAE_vector;
tspan = (0.0, 15.0);

##Find equilibrium point
inif! = (out,x) -> system_model!(out, dx0 ,x, (u0,OMIB2), 0.0)
sys_solve = nlsolve(inif!, x0)
x0_init = sys_solve.zero

##Construct DAE problem
prob = DiffEqBase.DAEProblem(system_model!, dx0, x0_init, tspan,
                            (u0, OMIB2), differential_vars = diff_vars)

##Solve and plot DAE
sol = solve(prob,IDA())
fig1 = plot(sol)

##Construct Callback
tstop = [1.0]
cb = DiffEqBase.DiscreteCallback(change_t_one, step_change3!)

##Solve problem with Callback (increase in Pd at t=1.0 from 0.4 to 0.5)
sol2 = solve(prob2, IDA(init_all = :false), callback=cb, tstops=tstop)
fig2 = plot(sol2)

##Save Figure
savefig(fig2, "plots/sims/OMIB_2states_wloads.png")

##Return turbine governor power to normal state
OMIB2.dyn_injections[1].P_ref = 0.4
