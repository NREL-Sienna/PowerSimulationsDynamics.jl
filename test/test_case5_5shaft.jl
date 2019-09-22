"""
Case 5:
This case study a three bus system with 1 machine located at bus 2.
The generator uses the model of a one d- one q- machine, and has a 5-mass shaft and a turbine governor.
The fault disconnects a circuit between buses 1 and 2, doubling its impedance.
"""

##################################################
############### LOAD DATA ########################
##################################################

############### Data Network ########################

nodes_case5   =    [Bus(1 , "Bus 1"  , "REF" , 0 , 1.05  , (min=0.94, max=1.06), 138),
                    Bus(2 , "Bus 2"  , "PV" , 0 , 1.00 , (min=0.94, max=1.06), 138),
                    Bus(3 , "Bus 3"  , "PV" , 0 , 1.00 , (min=0.94, max=1.06), 138)]

branch_case5 =     [Line("Line1", true, 0.0, 0.0, Arc(from=nodes_case5[1], to=nodes_case5[2]), 0.01, 0.05, (from=0.0, to=0.0), 100, 1.04),
                    Line("Line2", true, 0.0, 0.0, Arc(from=nodes_case5[2], to=nodes_case5[3]), 0.02, 0.1, (from=0.0, to=0.0), 100, 1.04),
                    Line("Line3", true, 0.0, 0.0, Arc(from=nodes_case5[1], to=nodes_case5[3]), 0.02, 0.1, (from=0.0, to=0.0), 100, 1.04)]

#Trip of a single circuit of Line 1 -> Resistance and Reactance doubled.
branch_case5_fault =     [Line("Line1", true, 0.0, 0.0, Arc(from=nodes_case5[1], to=nodes_case5[2]), 0.02, 0.1, (from=0.0, to=0.0), 100, 1.04),
                          Line("Line2", true, 0.0, 0.0, Arc(from=nodes_case5[2], to=nodes_case5[3]), 0.02, 0.1, (from=0.0, to=0.0), 100, 1.04),
                          Line("Line3", true, 0.0, 0.0, Arc(from=nodes_case5[1], to=nodes_case5[3]), 0.02, 0.1, (from=0.0, to=0.0), 100, 1.04)]

loads_case5 = [PowerLoad("Bus2", true, nodes_case5[2], PowerSystems.ConstantPower, 0.3, 0.05, 0.3, 0.05),
               PowerLoad("Bus3", true, nodes_case5[3], PowerSystems.ConstantPower, 0.3, 0.05, 0.3, 0.05)]


############### Data devices ########################

inf_gen_case5 = StaticSource(1, #number
                :InfBus, #name
                nodes_case5[1], #bus
                1.00, #VR
                0.0, #VI
                0.000001) #Xth

######## Machine Data #########

### Case 5: 4th Order Model with AVR + TG + Multishaft ###
case5_machine =  OneDOneQMachine(0.0, #R
                                  1.3125, #Xd
                                  1.2578, #Xq
                                  0.1813, #Xd_p
                                  0.25, #Xq_p
                                  5.89, #Td0_p
                                  0.6, #Tq0_p
                                  100.0)   #MVABase

######## Shaft Data #########

case5_shaft = FiveMassShaft(3.01, #5.148, #H
                            0.3348, #H_hp
                            0.7306, #H_ip
                            0.8154, #H_lp
                            0.0452, #H_ex,
                            0.0, #2.0, #D
                            0.5180, #D_hp
                            0.2240, #D_ip
                            0.2240, #D_lp
                            0.1450, #D_ex
                            0.0518, #D_12
                            0.0224, #D_23
                            0.0224, #D_34
                            0.0145, #D_45
                            33.07, #K_hp
                            28.59, #K_ip
                            44.68, #K_lp
                            21.984) #K_ex


######## PSS Data #########
cases_no_pss = PSSFixed(0.0)


######## TG Data #########

### TG for Case 5 ###
case5_tg = TGTypeII(0.05, #R
                    1.0, #T1
                    2.0, #T2
                    1.5, #τ_max
                    0.1) #τ_min

########  AVR Data #########
### AVRs for Case 2, 3, 4 and 5 ###
case5_avr =    AVRTypeI(20.0, #Ka - Gain
                        0.01, #Ke
                        0.063, #Kf
                        0.2, #Ta
                        0.314, #Te
                        0.35, #Tf
                        0.001, #Tr
                        5.0, #Vrmax
                        -5.0, #Vrmin
                        0.0039, #Ae - 1st ceiling coefficient
                        1.555) #Be - 2nd ceiling coefficient



### Case 5 Generator ###
case5_gen = DynGenerator(1, #Number
                         :Case5Gen,
                         nodes_case5[2], #bus
                         1.0, # ω_ref,
                         1.0155, #V_ref
                         0.5, #P_ref
                         case5_machine, #machine
                         case5_shaft, #shaft
                         case5_avr, #avr
                         case5_tg, #tg
                         cases_no_pss) #pss


######################### Dynamical System ########################

case5_DynSystem = DynamicSystem(nodes_case5,
                                branch_case5,
                                [case5_gen],
                                vcat(inf_gen_case5,loads_case5),
                                100.0,
                                60.0)


##################################################
############### SOLVE PROBLEM ####################
##################################################


#Compute Y_bus after fault
Ybus_fault = PSY.Ybus(branch_case5_fault, nodes_case5)[:,:];

#Initialize variables
dx0 = zeros(LITS.get_total_rows(case5_DynSystem))
x0 = [1.02, 1.0, 1.0, 0.0, -0.01, -0.01,
      1.0, #eq_p
      0.0, #ed_p
      0.05, #δ
      1.0, #ω
      0.05, #δ_HP
      1.0, #ω_HP
      0.05, #δ_IP
      1.0, #ω_LP
      0.05, #δ_LP
      1.0, #ω_LP
      0.05, #δ_Ex
      1.0, #ω_Ex
      1.0, #Vf
      0.01, #Vr1
      -0.1, #Vr2,
      1.0, #Vm
      0.0] #xg
diff_vars = case5_DynSystem.DAE_vector
tspan = (0.0, 20.0);

#Find initial condition
inif! = (out,x) -> system_model!(out, dx0 ,x, (Ybus_fault,case5_DynSystem), 0.0)
sys_solve = nlsolve(inif!, x0)
x0_init = sys_solve.zero

#Define problem
prob = DiffEqBase.DAEProblem(system_model!, dx0, x0_init, tspan,
                            (Ybus_fault, case5_DynSystem), differential_vars = diff_vars);

#Solve problem in equilibrium
sol = solve(prob, IDA());

#Define data for using callbacks for defining the fault
tstop = [1.0] #Define a timestop at t=1, the step change
cb = DiffEqBase.DiscreteCallback(LITS.change_t_one, LITS.Y_change!)

#Solve DAE system
sol2 = solve(prob, IDA(init_all = :false), dtmax= 0.02, callback=cb, tstops=tstop);

#Obtain data for angles
series = LITS.get_state_series(sol2, case5_DynSystem, (:Case5Gen, :δ));

@test sol2.retcode == :Success
