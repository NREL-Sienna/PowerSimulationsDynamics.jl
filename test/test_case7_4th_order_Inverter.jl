"""
Case 7:
This case study a three bus system with 1 machine (One d- One q-: 4th order model), a VSM of 19 states and an infinite source.
The perturbation increase the reference power (analogy for mechanical power) of the machine from 0.6 to 0.8.
"""


##################################################
############### LOAD DATA ########################
##################################################

############### Data Network ########################

nodes_case7   = [Bus(1 , "Bus 1"  , "REF" , 0 , 1.02  , (min=0.94, max=1.06), 138),
                    Bus(2 , "Bus 2"  , "PV" , 0 , 1.00 , (min=0.94, max=1.06), 138),
                    Bus(3 , "Bus 3"  , "PQ" , 0 , 1.00 , (min=0.94, max=1.06), 138)]

branch_case7  =  [Line("Line1", true, 0.0, 0.0, Arc(from=nodes_case7[1], to=nodes_case7[3]), 0.01, 0.12, (from=0.0, to=0.0), 100, 1.04),
                    Line("Line2", true, 0.0, 0.0, Arc(from=nodes_case7[1], to=nodes_case7[2]), 0.01, 0.12, (from=0.0, to=0.0), 100, 1.04),
                    Line("Line3", true, 0.0, 0.0, Arc(from=nodes_case7[2], to=nodes_case7[3]), 0.01, 0.12, (from=0.0, to=0.0), 100, 1.04)]

#Trip of Line 1.
branch_case7_fault = [Line("Line2", true, 0.0, 0.0, Arc(from=nodes_case7[1], to=nodes_case7[2]), 0.01, 0.12, (from=0.0, to=0.0), 100, 1.04),
                        Line("Line3", true, 0.0, 0.0, Arc(from=nodes_case7[2], to=nodes_case7[3]), 0.01, 0.12, (from=0.0, to=0.0), 100, 1.04)]


loads_case7 =   [PowerLoad("Bus1", true, nodes_case7[1], PowerSystems.ConstantPower, 0.5, 0.1, 1.5, 0.8),
                   PowerLoad("Bus2", true, nodes_case7[2], PowerSystems.ConstantPower, 1.0, 0.3, 1.5, 0.8),
                   PowerLoad("Bus3", true, nodes_case7[3], PowerSystems.ConstantPower, 0.3, 0.1, 0.5, 0.3)]


############### Data devices ########################

inf_gen_case7 = StaticSource(1, #number
                    :InfBus, #name
                    nodes_case7[1],#bus
                    1.00, #VR
                    0.0, #VI
                    0.000005) #Xth


######## Machine Data #########

### Case 2: 4th Order Model with AVR (3-bus case) ###
case7_machine =  OneDOneQMachine(0.0, #R
                                  1.3125, #Xd
                                  1.2578, #Xq
                                  0.1813, #Xd_p
                                  0.25, #Xq_p
                                  5.89, #Td0_p
                                  0.6, #Tq0_p
                                  100.0)   #MVABase

######## Shaft Data #########

### Shafts for Gen ###
case7_shaft = SingleMass(3.01, #H (M = 6.02 -> H = M/2)
                            0.0) #D

######## PSS Data #########
cases_no_pss = PSSFixed(0.0)


######## TG Data #########

### No TG for Cases 1, 2, 3, 4 ###
case7_no_tg = TGFixed(1.0) #eff


########  AVR Data #########
### AVRs for Case 2, 3, 4 and 5 ###
case7_avr = AVRTypeI(20.0, #Ka - Gain
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

### Case 7 Generators ###
case7_gen = DynGenerator(1, #Number
                         :Case7Gen,
                         nodes_case7[2], #bus
                         1.0, # ω_ref,
                         1.0142, #V_ref
                         0.6, #P_ref
                         case7_machine, #machine
                         case7_shaft, #shaft
                         case7_avr, #avr
                         case7_no_tg, #tg
                         cases_no_pss) #pss

############### Inverter Data ########################

converter = AvgCnvFixedDC(138.0, #Rated Voltage
                          100.0) #Rated MVA

dc_source = FixedDCSource(1500.0) #Not in the original data, guessed.

filt = LCLFilter(0.08, #Series inductance lf in pu
                   0.003, #Series resitance rf in pu
                   0.074, #Shunt capacitance cf in pu
                   0.2, #Series reactance rg to grid connection (#Step up transformer or similar)
                   0.01) #Series resistance lg to grid connection (#Step up transformer or similar)

pll = PLL(500.0, #ω_lp: Cut-off frequency for LowPass filter of PLL filter.
          0.084, #k_p: PLL proportional gain
          4.69) #k_i: PLL integral gain

virtual_H = VirtualInertia(2.0, #Ta:: VSM inertia constant
                           400.0, #kd:: VSM damping coefficient
                           20.0, #kω:: Frequency droop gain in pu
                           2*pi*50.0) #ωb:: Rated angular frequency

Q_control = ReactivePowerDroop(0.2, #kq:: Reactive power droop gain in pu
                              1000.0) #ωf:: Reactive power cut-off low pass filter frequency

outer_control = VirtualInertiaQdroop(virtual_H, Q_control)

vsc = CombinedVIwithVZ(0.59, #kpv:: Voltage controller proportional gain
                       736.0, #kiv:: Voltage controller integral gain
                       0.0, #kffv:: Binary variable enabling the voltage feed-forward in output of current controllers
                       0.0, #rv:: Virtual resistance in pu
                       0.2, #lv: Virtual inductance in pu
                       1.27, #kpc:: Current controller proportional gain
                       14.3, #kiv:: Current controller integral gain
                       0.0, #kffi:: Binary variable enabling the current feed-forward in output of current controllers
                       50.0, #ωad:: Active damping low pass filter cut-off frequency
                       0.2) #kad:: Active damping gain

case7_inv = DynInverter(2, #number
                             :DARCO, #name
                             nodes_case7[3], #bus location
                             1.0, #ω_ref
                             1.02, #V_ref
                             0.5, #P_ref
                             0.0, #Q_ref
                             100.0, #MVABase
                             converter, #Converter
                             outer_control, #OuterControl
                             vsc, #Voltage Source Controller
                             dc_source, #DC Source
                             pll, #Frequency Estimator
                             filt) #Output Filter


######################### Dynamical System ########################

case7_DynSystem = DynamicSystem(nodes_case7, branch_case7, [case7_inv, case7_gen], vcat(inf_gen_case7, loads_case7), 100.0, 50.0);


##################################################
############### SOLVE PROBLEM ####################
##################################################

dx0 = zeros(LITS.get_total_rows(case7_DynSystem))
x0 = [1.00, #V1_R
          1.00, #V2_R
          1.00, #V3_R
          0.0, #V1_I
         -0.01, #V2_I
         -0.01, #V3_I
         0.0, #δω_vsm
          0.2, #δθ_vsm
          0.025, #qm
          0.0015, #ξ_d
          -0.07, #ξ_q
          0.05, #γ_d
        -0.001, #γ_q
         0.95, #ϕ_d
         -0.10, #ϕ_q
         1.004, #vpll_d
         0.0, #vpll_q
         0.0, #ε_pll
         0.1, #δθ_pll
         0.5, #id_cv
         0.0, #iq_cv
         0.95, #vod
         -0.1, #voq
         0.49, #iod
        -0.1, #ioq
          1.0, #eq_p
          0.47, #ed_p
          0.6, #δ
          1.0, #ω
          2.1, #Vf
          0.28, #Vr1
          -0.39, #Vr2,
          1.0] #Vm

diff_vars = case7_DynSystem.DAE_vector
case7_inv.inner_vars[13] = 0.95 #Vd_cnv var
case7_inv.inner_vars[14] = -0.1 #Vq_cnv var
Ybus_fault = Ybus(branch_case7_fault, nodes_case7)[:,:]
u0 = [0.2]
tspan = (0.0, 20.0);


#Find initial condition
inif! = (out,x) -> system_model!(out, dx0 ,x, (Ybus_fault,case7_DynSystem), 0.0)
sys_solve = nlsolve(inif!, x0, xtol=:1e-8,ftol=:1e-8,method=:trust_region)
x0_init = sys_solve.zero

#Define problem
prob = DiffEqBase.DAEProblem(system_model!, dx0, x0_init, tspan, (u0, case7_DynSystem), differential_vars = diff_vars)

#Solve problem in equilibrium
sol = solve(prob, IDA());

#Define data for using callbacks for defining the fault
tstop = [1.0] #Define a timestop at t=1, the step change
cb = DiffEqBase.DiscreteCallback(LITS.change_t_one, LITS.step_change4!)

#Solve DAE system
sol2 = solve(prob, IDA(), dtmax= 0.02, callback=cb, tstops=tstop)

#Obtain data for angles
series = LITS.get_state_series(sol2, case1_DynSystem, (:Case1Gen, :δ))

#Obtain voltage magnitudes at inverter bus
series2 = LITS.get_voltagemag_series(sol2, case7_DynSystem, 3)

@test sol2.retcode == :Success
