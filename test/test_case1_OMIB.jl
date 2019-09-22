"""
Case 1:
This case study defines a classical machine against an infinite bus. The fault
drop a circuit on the (double circuit) line connecting the two buses, doubling its impedance
"""


##################################################
############### LOAD DATA ########################
##################################################

############### Data Network ########################

nodes_case1 = [Bus(1 , #number
                   "Bus 1", #Name
                   "REF" , #BusType (REF, PV, PQ)
                   0, #Angle in radians
                   1.05, #Voltage in pu
                   (min=0.94, max=1.06), #Voltage limits in pu
                   69), #Base voltage in kV
                   Bus(2 , "Bus 2"  , "PV" , 0 , 1.0 , (min=0.94, max=1.06), 69)]

branch_case1 = [Line("Line1", #name
                     true, #available
                     0.0, #active power flow initial condition (from-to)
                     0.0, #reactive power flow initial condition (from-to)
                     Arc(from=nodes_case1[1], to=nodes_case1[2]), #Connection between buses
                     0.01, #resistance in pu
                     0.05, #reactance in pu
                     (from=0.0, to=0.0), #susceptance in pu
                     18.046, #rate in MW
                     1.04)]  #angle limits (-min and max)

#Trip of a single circuit of Line 1 -> Resistance and Reactance doubled.
branch_case1_fault = [Line("Line1", #name
                           true, #available
                           0.0, #active power flow initial condition (from-to)
                           0.0, #reactive power flow initial condition (from-to)
                           Arc(from=nodes_case1[1], to=nodes_case1[2]), #Connection between buses
                           0.02, #resistance in pu
                           0.1, #reactance in pu
                           (from=0.0, to=0.0), #susceptance in pu
                           18.046, #rate in MW
                           1.04)]  #angle limits (-min and max)

loads_case1 = [PowerLoad("Bus1", true, nodes_case1[2], PowerSystems.ConstantPower, 0.3, 0.01, 0.3, 0.01)]


############### Data devices ########################

inf_gen_case1 = StaticSource(1, #number
                :InfBus, #name
                nodes_case1[1], #bus
                1.05, #VR
                0.0, #VI
                0.000001) #Xth

######## Machine Data #########

### Case 1: Classical machine against infinite bus ###
case1_machine = BaseMachine(0.0, #R
                            0.2995, #Xd_p
                            0.7087, #eq_p
                            100.0)  #MVABase

######## Shaft Data #########

### Shaft for Case 1 ###
case1_shaft = SingleMass(3.148, #H
                     2.0) #D



######## PSS Data #########
cases_no_pss = PSSFixed(0.0)


######## TG Data #########

### No TG for Cases 1, 2, 3, 4 ###
case1234_no_tg = TGFixed(1.0) #eff

########  AVR Data #########
case1_avr = AVRFixed(0.0) #Vf not applicable in Classic Machines

### Case 1 Generator ###
case1_gen = DynGenerator(1, #Number
                      :Case1Gen,
                      nodes_case1[2], #bus
                      1.0, # ω_ref,
                      1.0, #V_ref
                      0.5, #P_ref
                      case1_machine, #machine
                      case1_shaft, #shaft
                      case1_avr, #avr
                      case1234_no_tg, #tg
                      cases_no_pss) #pss


######################### Dynamical System ########################

case1_DynSystem = DynamicSystem(nodes_case1,
                              branch_case1,
                              [case1_gen],
                              vcat(inf_gen_case1,loads_case1),
                              100.0,
                              60.0)


##################################################
############### SOLVE PROBLEM ####################
##################################################

#Compute Y_bus after fault
Ybus_fault = PSY.Ybus(branch_case1_fault, nodes_case1)[:,:]

#Initialize variables
dx0 = zeros(LITS.get_total_rows(case1_DynSystem))
x0 = [1.05, #VR_1
      1.0, #VR_21.0, 0.0, -0.01, -0.01,
      0.0, #VI_1
      0.01, #VI_2
      0.2, #δ
      1.0] #ω
diff_vars = case1_DynSystem.DAE_vector
tspan = (0.0, 20.0);

#Find initial condition
inif! = (out,x) -> system_model!(out, dx0 ,x, (Ybus_fault,case1_DynSystem), 0.0)
sys_solve = nlsolve(inif!, x0)
x0_init = sys_solve.zero

#Define problem
prob = DiffEqBase.DAEProblem(system_model!, dx0, x0_init, tspan,
                            (Ybus_fault, case1_DynSystem), differential_vars = diff_vars)

#Solve problem in equilibrium
sol = solve(prob, IDA());

#Define data for using callbacks for defining the fault
tstop = [1.0] #Define a timestop at t=1, the step change
cb = DiffEqBase.DiscreteCallback(LITS.change_t_one, LITS.Y_change!)

#Solve DAE system
sol2 = solve(prob, IDA(init_all = :false), dtmax= 0.02, callback=cb, tstops=tstop)

#Obtain data for angles
series = LITS.get_state_series(sol2, case1_DynSystem, (:Case1Gen, :δ))

@test sol2.retcode == :Success
