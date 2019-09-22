"""
Sources data
"""

inf_gen_case1 = StaticSource(1, #number
                :InfBus, #name
                nodes_case1[1], #bus
                1.05, #VR
                0.0, #VI
                0.000001) #Xth

inf_gen_case234 = StaticSource(1, #number
                :InfBus, #name
                nodes_case234[1], #bus
                1.02, #VR
                0.0, #VI
                0.000001) #Xth

inf_gen_case5 = StaticSource(1, #number
                :InfBus, #name
                nodes_case5[1], #bus
                1.00, #VR
                0.0, #VI
                0.000001) #Xth


"""
Generators data
"""

######## Machine Data #########

### Case 1: Classical machine against infinite bus ###
case1_machine = BaseMachine(0.0, #R
                            0.2995, #Xd_p
                            0.7087, #eq_p
                            100.0)  #MVABase



### Case 2: 4th Order Model with AVR (3-bus case) ###
case2_machine2 =  OneDOneQMachine(0.0, #R
                                  1.3125, #Xd
                                  1.2578, #Xq
                                  0.1813, #Xd_p
                                  0.25, #Xq_p
                                  5.89, #Td0_p
                                  0.6, #Tq0_p
                                  100.0)   #MVABase

case2_machine3 =  OneDOneQMachine(0.0, #R
                                  1.3125, #Xd
                                  1.2578, #Xq
                                  0.1813, #Xd_p
                                  0.25, #Xq_p
                                  5.89, #Td0_p
                                  0.6, #Tq0_p
                                  100.0)   #MVABase

### Case 3: 6th Order Model with AVR (3-bus case) ###
case3_machine2 = SimpleMarconatoMachine(0.0,
                                        1.3125, #Xd
                                        1.2578, #Xq
                                        0.1813, #Xd_p
                                        0.25, #Xq_p
                                        0.14, #Xd_pp
                                        0.18, #Xq_pp
                                        5.89, #Td0_p
                                        0.6, #Tq0_p
                                        0.5, #Td0_pp
                                        0.023, #Tq0_pp
                                        0.0, #T_AA
                                        100.0) #MVABase

case3_machine3 = SimpleMarconatoMachine(0.0,
                                        1.3125, #Xd
                                        1.2578, #Xq
                                        0.1813, #Xd_p
                                        0.25, #Xq_p
                                        0.14, #Xd_pp
                                        0.18, #Xq_pp
                                        5.89, #Td0_p
                                        0.6, #Tq0_p
                                        0.5, #Td0_pp
                                        0.023, #Tq0_pp
                                        0.0, #T_AA
                                        100.0) #MVABase

### Case 4: 8th Order Model with AVR (3-bus case) ###
case4_machine2 = MarconatoMachine(0.0,
                                  1.3125, #Xd
                                  1.2578, #Xq
                                  0.1813, #Xd_p
                                  0.25, #Xq_p
                                  0.14, #Xd_pp
                                  0.18, #Xq_pp
                                  5.89, #Td0_p
                                  0.6, #Tq0_p
                                  0.5, #Td0_pp
                                  0.023, #Tq0_pp
                                  0.0, #T_AA
                                  100.0) #MVABase


case4_machine3 = MarconatoMachine(0.0,
                                  1.3125, #Xd
                                  1.2578, #Xq
                                  0.1813, #Xd_p
                                  0.25, #Xq_p
                                  0.14, #Xd_pp
                                  0.18, #Xq_pp
                                  5.89, #Td0_p
                                  0.6, #Tq0_p
                                  0.5, #Td0_pp
                                  0.023, #Tq0_pp
                                  0.0, #T_AA
                                  100.0) #MVABase

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

### Shaft for Case 1 ###
case1_shaft = SingleMass(3.148, #H
                         2.0) #D

### Shafts for Gens 2 and 3: Cases 2, 3 and 4 ###
case234_shaft2 = SingleMass(3.01, #H (M = 6.02 -> H = M/2)
                            0.0) #D

case234_shaft3 = SingleMass(3.01, #H (M = 6.02 -> H = M/2)
                            0.0) #D

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

### No TG for Cases 1, 2, 3, 4 ###
case1234_no_tg = TGFixed(1.0) #eff

### TG for Case 5 ###
case5_tg = TGTypeII(0.05, #R
                    1.0, #T1
                    2.0, #T2
                    1.5, #τ_max
                    0.1) #τ_min

######## AVR Data #########

### AVR for Case 1 ###
case1_avr = AVRFixed(0.0) #Vf not applicable in Classic Machines

### AVRs for Case 2, 3, 4 and 5 ###
case2345_avr2 = AVRTypeI(20.0, #Ka - Gain
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

case2345_avr3 = AVRTypeI(20.0, #Ka - Gain
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





######################### Generators ########################

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


### Case 2 Generators ###
case2_gen2 = DynGenerator(1, #Number
                         :Case2Gen2,
                         nodes_case234[2], #bus
                         1.0, # ω_ref,
                         1.0142, #V_ref
                         1.0, #P_ref
                         case2_machine2, #machine
                         case234_shaft2, #shaft
                         case2345_avr2, #avr
                         case1234_no_tg, #tg
                         cases_no_pss) #pss

case2_gen3 = DynGenerator(2, #Number
                         :Case2Gen3,
                         nodes_case234[3], #bus
                         1.0, # ω_ref,
                         1.0059, #V_ref
                         1.0, #P_ref
                         case2_machine3, #machine
                         case234_shaft3, #shaft
                         case2345_avr3, #avr
                         case1234_no_tg, #tg
                         cases_no_pss) #pss

### Case 3 Generators ###
case3_gen2 = DynGenerator(1, #Number
                         :Case3Gen2,
                         nodes_case234[2], #bus
                         1.0, # ω_ref,
                         1.0142, #V_ref
                         1.0, #P_ref
                         case3_machine2, #machine
                         case234_shaft2, #shaft
                         case2345_avr2, #avr
                         case1234_no_tg, #tg
                         cases_no_pss) #pss

case3_gen3 = DynGenerator(2, #Number
                         :Case3Gen3,
                         nodes_case234[3], #bus
                         1.0, # ω_ref,
                         1.0059, #V_ref
                         1.0, #P_ref
                         case3_machine3, #machine
                         case234_shaft3, #shaft
                         case2345_avr3, #avr
                         case1234_no_tg, #tg
                         cases_no_pss) #pss


### Case 4 Generators ###
case4_gen2 = DynGenerator(1, #Number
                         :Case4Gen2,
                         nodes_case234[2], #bus
                         1.0, # ω_ref,
                         1.0142, #V_ref
                         1.0, #P_ref
                         case4_machine2, #machine
                         case234_shaft2, #shaft
                         case2345_avr2, #avr
                         case1234_no_tg, #tg
                         cases_no_pss) #pss

case4_gen3 = DynGenerator(2, #Number
                         :Case4Gen3,
                         nodes_case234[3], #bus
                         1.0, # ω_ref,
                         1.0059, #V_ref
                         1.0, #P_ref
                         case4_machine3, #machine
                         case234_shaft3, #shaft
                         case2345_avr3, #avr
                         case1234_no_tg, #tg
                         cases_no_pss) #pss

### Case 5 Generator ###
case5_gen = DynGenerator(1, #Number
                         :Case5Gen,
                         nodes_case5[2], #bus
                         1.0, # ω_ref,
                         1.0155, #V_ref
                         0.5, #P_ref
                         case5_machine, #machine
                         case5_shaft, #shaft
                         case2345_avr2, #avr
                         case5_tg, #tg
                         #case1234_no_tg,
                         cases_no_pss) #pss




######################### Dynamical System ########################

case1_DynSystem = DynamicSystem(nodes_case1,
                                branch_case1,
                                [case1_gen],
                                vcat(inf_gen_case1,loads_case1),
                                100.0,
                                60.0)


case2_DynSystem = DynamicSystem(nodes_case234,
                                branch_case234,
                                [case2_gen2, case2_gen3],
                                vcat(inf_gen_case234,loads_case234),
                                100.0,
                                60.0)


case3_DynSystem = DynamicSystem(nodes_case234,
                                branch_case234,
                                [case3_gen2, case3_gen3],
                                vcat(inf_gen_case234,loads_case234),
                                100.0,
                                60.0)

case4_DynSystem = DynamicSystem(nodes_case234,
                                branch_case234,
                                [case4_gen2, case4_gen3],
                                vcat(inf_gen_case234,loads_case234),
                                100.0,
                                60.0)

case5_DynSystem = DynamicSystem(nodes_case5,
                                branch_case5,
                                [case5_gen],
                                vcat(inf_gen_case5,loads_case5),
                                100.0,
                                60.0)
