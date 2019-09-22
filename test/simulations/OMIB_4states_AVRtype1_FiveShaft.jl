##Construct Dynamical System
OMIB7 = DynamicSystem(nodes_OMIB, #2bus Network
                      branch_OMIB, #Network data
                      [Gen6AVR], #4StateGen with AVR and 5 Shaft
                      [inf_gen], #Inf bus and load
                      100.0, #MVABase
                      60.0) #f_sys

##Initial conditions
dx0 = zeros(OMIB7.counts[:total_rows])
x0 = [1.1, #v_1R
    1.0, #v_2R
    0.0, #v_1I
    0.0, #v_2I
    1.05, #eq_p
    0.10, #ed_pp
    0.0, #δ
    1.0, #ω
    0.0, #δ_hp
    1.0, #ω_hp
    0.0, #δ_ip
    1.0, #ω_ip
    0.0, #δ_lp
    1.0, #ω_lp
    0.0, #δ_ex
    1.0, #ω_ex
    1.05, #Vf AVR
    0.5, #Vr1 AVR
    0.0, #Vr2 AVR
    1.05] #Vm AVR
u0 = [0.4] #Pd
diff_vars = OMIB7.DAE_vector;
tspan = (0.0, 20.0)

##Find equilibrium point
inif! = (out,x) -> system_model!(out, dx0 ,x, (u0,OMIB7), 0.0)
sys_solve7 = nlsolve(inif!, x0)
x0_init = sys_solve7.zero

##Construct DAE problem
prob = DiffEqBase.DAEProblem(system_model!, dx0, x0_init, tspan,
                            (u0, OMIB7), differential_vars = diff_vars)

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
savefig(fig2, "plots/sims/OMIB_4states_AVRtype1_FiveShaft.png")

##Return turbine governor power to normal state
OMIB7.dyn_injections[1].P_ref = 0.4
