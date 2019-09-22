# Electrical Parameters of Machine 1, Milano's Book Page 526
################## Machine Data #####################
Basic = BaseMachine(   0.0, #R
                       0.2995, #Xd_p
                       1.05, #eq_p
                       615.0)  #MVABase

oneDoneQ = OneDOneQMachine(0.0, #R
                        0.8979, #Xd
                        0.646, #Xq
                        0.2995, #Xd_p
                        0.04, #Xq_p
                        7.4, #Td0_p
                        0.033, #Tq0_p
                        615.0)   #MVABase

AndersonFouad = AndersonFouadMachine(0.0, #R
                        0.8979, #Xd
                        0.646, #Xq
                        0.2995, #Xd_p
                        0.646, #Xq_p
                        0.23, #Xd_pp
                        0.4, #Xq_pp
                        7.4, #Td0_p
                        0.01, #Tq0_p #Data not available in Milano: Used 0.01
                        0.03, #Td0_pp
                        0.033, #Tq0_pp
                        615.0) #MVABase

KundurMachine = SimpleFullMachine(0.003, #R on Example 3.1 and 4.1 of Kundur
                                  0.0006, #R_f
                                  0.0284, #R_1d or RD in Machowski
                                  0.0062, #R_1q or RQ on Machowski
                                  1.81, #L_d
                                  1.76, #L_q
                                  1.66, #L_ad or k*M_f or k*M_D in Machowski
                                  1.61, #L_aq or k*M_Q in Machowski
                                  1.66, #L_f1d or L_fD in Machowski. Assumed to be equal to L_ad
                                  1.825, #L_ff
                                  0.1713, #L_1d or L_D in Machowski
                                  0.7525, #L_1q or L_Q in Machowski
                                  555.0) #MVABase

KundurFullMachine = FullMachine(0.003, #R on Example 3.1 and 4.1 of Kundur
                                0.0006, #R_f
                                #0.003, #R_f
                                0.0284, #R_1d or RD in Machowski
                                0.0062, #R_1q or RQ on Machowski
                                1.81, #L_d
                                1.76, #L_q
                                1.66, #L_ad or k*M_f or k*M_D in Machowski
                                1.61, #L_aq or k*M_Q in Machowski
                                1.66, #L_f1d or L_fD in Machowski. Assumed to be equal to L_ad
                                1.825, #L_ff
                                0.1713, #L_1d or L_D in Machowski
                                0.7525, #L_1q or L_Q in Machowski
                                555.0) #MVABase

Mach2_benchmark = OneDOneQMachine(0.0, #R
                        1.3125, #Xd
                        1.2578, #Xq
                        0.1813, #Xd_p
                        0.25, #Xq_p
                        5.89, #Td0_p
                        0.6, #Tq0_p
                        100.0)   #MVABase

Mach3_benchmark = OneDOneQMachine(0.0, #R
                        1.3125, #Xd
                        1.2578, #Xq
                        0.1813, #Xd_p
                        0.25, #Xq_p
                        5.89, #Td0_p
                        0.6, #Tq0_p
                        100.0)   #MVABase


################ Shaft Data #####################
BaseShaft = SingleMass(5.148, #H
                       2.0) #D

FiveShaft = FiveMassShaft(5.148,  #H
                          0.3348, #H_hp
                          0.7306, #H_ip
                          0.8154, #H_lp
                          0.0452, #H_ex,
                          2.0,    #D
                          0.5180, #D_hp
                          0.2240, #D_ip
                          0.2240, #D_lp
                          0.1450, #D_ex
                          0.0518, #D_12
                          0.0224, #D_23
                          0.0224, #D_34
                          0.0145, #D_45
                          33.07,  #K_hp
                          28.59,  #K_ip
                          44.68,  #K_lp
                          21.984) #K_ex

KundurShaft = SingleMass(3.525, #H given in example 3.5 in Kundur
                         2.0) #D: Damping is not given, so I assume the value 2.0

Shaft2_benchmark = SingleMass(3.01, #H (M = 6.02 -> H = M/2)
                       0.0) #D

Shaft3_benchmark = SingleMass(3.01, #H (M = 6.02 -> H = M/2)
                      0.0) #D

################# PSS Data #####################
no_pss = PSSFixed(0.0)

################ TG Data #####################
fixed_tg = TGFixed(1.0) #eff

typeI_tg = TGTypeI(0.02, #R
                   0.1, #Ts
                   0.45, #Tc
                   0.0, #T3
                   0.0, #T4
                   50.0, #T5
                   0.3, #P_min
                   1.2) #P_max

typeII_tg = TGTypeII(0.05, #R
                     0.3, #T1
                     0.1, #T2
                     1.0, #τ_max
                     0.1) #τ_min

################ AVR Data #####################
proportional_avr = AVRSimple(5000.0) #Kv

fixed_avr = AVRFixed(1.05) #Emf

typeI_avr = AVRTypeI(200.0, #Ka
                     1.0, #Ke
                     0.0012, #Kf
                     0.02, #Ta
                     0.19, #Te
                     1.0, #Tf
                     0.001, #Tr
                     9.9, #Vr_max
                     0.0, #Vr_min
                     0.0006, #Ae
                     0.9)

kundur_avr = AVRFixed(0.00043) #Emf
gen_avr_bench = AVRFixed(1.7613) #Vf

gen2_avr_benchmark_m = AVRTypeI(20.0, #Ka - Gain
                   0.01, #Ke
                   0.063, #Kf
                   0.2, #Ta
                   0.314, #Te
                   0.35, #Tf
                   0.001, #Tr
gen2_avr_benchmark = AVRTypeII(20.0, #K0 - Gain
                     0.2, #T1 - 1st pole
                     0.063, #T2 - 1st zero
                     0.35, #T3 - 2nd pole
                     0.01, #T4 - 2nd zero
                     0.314, #Te - Field current time constant
                     0.001, #Tr - Measurement time constant
                     5.0, #Vrmax
                     -5.0, #Vrmin
                     0.0039, #Ae - 1st ceiling coefficient
                     1.555) #Be - 2nd ceiling coefficient

gen3_avr_benchmark = AVRTypeII(20.0, #K0 - Gain
                  0.2, #T1 - 1st pole
                  0.063, #T2 - 1st zero
                  0.35, #T3 - 2nd pole
                  0.01, #T4 - 2nd zero
                  0.314, #Te - Field current time constant
                  0.001, #Tr - Measurement time constant
                  5.0, #Vrmax
                  -5.0, #Vrmin
                  0.0039, #Ae - 1st ceiling coefficient
                  1.555) #Be - 2nd ceiling coefficient

gen2_avr_benchmark_m = AVRTypeIIManual(20.0, #Ka - Gain
                   0.2, #Ta - 1st pole
                   0.063, #Kf - 1st zero
                   0.35, #Tf - 2nd pole
                   0.01, #Te - 2nd zero
                   0.314, #Tr - Field current time constant
                   5.0, #Vrmax
                   -5.0, #Vrmin
                   0.0039, #Ae - 1st ceiling coefficient
                   1.555) #Be - 2nd ceiling coefficient
gen3_avr_benchmark_m = AVRTypeI(20.0, #Ka - Gain
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

gen2_avr_fixed = AVRFixed(2.1996)
gen3_avr_fixed = AVRFixed(1.7327)
gen3_avr_benchmark_m = AVRTypeIIManual(20.0, #Ka - Gain
                0.2, #Ta - 1st pole
                0.063, #Kf - 1st zero
                0.35, #Tf - 2nd pole
                0.01, #Te - 2nd zero
                0.314, #Tr - Field current time constant
                5.0, #Vrmax
                -5.0, #Vrmin
                0.0039, #Ae - 1st ceiling coefficient
                1.555) #Be - 2nd ceiling coefficient

######################### Generators ########################

Gen1AVR = DynGenerator(1, #Number
                 :TestGen,
                 nodes_OMIB[2],#bus
                 1.0, # ω_ref,
                 1.05,
                 0.4,
                 Basic,
                 BaseShaft,
                 proportional_avr, #avr
                 fixed_tg, #tg
                 no_pss)

Gen1AVRnoAVR = DynGenerator(1, #Number
                 :TestGen,
                 nodes_OMIB[2],#bus
                 1.0, # ω_ref,
                 1.05,
                 0.4,
                 Basic,
                 BaseShaft,
                 fixed_avr, #avr
                 fixed_tg, #tg
                 no_pss)

Gen2AVRnoAVR = DynGenerator(1, #Number
                 :TestGen,
                 nodes_OMIB[2],#bus
                 1.0, # ω_ref,
                 1.02,
                 0.4,
                 oneDoneQ,
                 BaseShaft,
                 fixed_avr, #avr
                 fixed_tg, #tg
                 no_pss)

Gen2AVR = DynGenerator(1, #Number
                 :TestGen,
                 nodes_OMIB[2],#bus
                 1.0, # ω_ref,
                 1.02,
                 0.4,
                 oneDoneQ,
                 BaseShaft,
                 proportional_avr, #avr
                 fixed_tg, #tg
                 no_pss)

Gen3AVR = DynGenerator(1, #number
                       :TestGen, #name
                       nodes_OMIB[2], #bus
                       1.0, # ω_ref
                       1.02,
                       0.4,
                       oneDoneQ, #machine
                       BaseShaft,
                       typeI_avr, #avr
                       fixed_tg, #tg
                       no_pss)

Gen4AVRTG = DynGenerator(1, #number
                       :TestGen, #name
                       nodes_OMIB[2], #bus
                       1.0, # ω_ref
                       1.02,
                       0.4,
                       oneDoneQ, #machine
                       BaseShaft,
                       typeI_avr, #avr
                       typeI_tg, #tg
                       no_pss)

Gen5AVR = DynGenerator(1, #number
                       :TestGen, #name
                       nodes_OMIB[2], #bus
                       1.0, # ω_ref
                       1.02,
                       0.4,
                       AndersonFouad, #machine
                       BaseShaft,
                       typeI_avr, #avr
                       fixed_tg, #tg
                       no_pss)

Gen6AVR = DynGenerator(1, #number
                       :TestGen, #name
                       nodes_OMIB[2], #bus
                       1.0, #ω_ref
                       1.02, #V_ref
                       0.4, #P_ref
                       oneDoneQ, #machine
                       FiveShaft, #shaft
                       typeI_avr, #avr
                       fixed_tg, #tg
                       no_pss)

GenKundur = DynGenerator(1, #number
                         :KundurGen, #name
                         nodes_kundur[1],
                         1.0, #ω_ref
                         1.02, #V_ref
                         0.2, #P_ref
                         KundurMachine, #machine
                         KundurShaft, #shaft
                         kundur_avr, #avr
                         fixed_tg, #tg
                         no_pss) #pss

GenFullKundur =  DynGenerator(1, #number
                              :KundurGen, #name
                              nodes_kundur[1],
                              1.0, #ω_ref
                              1.02, #V_ref
                              0.20, #P_ref
                              KundurFullMachine, #machine
                              KundurShaft, #shaft
                              kundur_avr, #avr
                              fixed_tg, #tg
                              no_pss) #pss


Gen2_benchmark = DynGenerator(1, #number
                              :Gen1, #name
                              nodes_benchmark[2],
                              1.0, #ω_ref
                              1.0142, #V_ref
                              1.0, #P_ref
                              Mach2_benchmark, #machine
                              Shaft2_benchmark, #shaft
                              gen2_avr_benchmark_m, #avr
                              fixed_tg, #tg
                              no_pss) #pss
Gen2_benchmark_tg = DynGenerator(1, #number
                              :Gen1, #name
                              nodes_benchmark[2],
                              1.0, #ω_ref
                              1.0142, #V_ref
                              1.0, #P_ref
                              Mach2_benchmark, #machine
                              Shaft2_benchmark, #shaft
                              gen2_avr_benchmark_m, #avr
                              typeII_tg, #tg
                              no_pss) #pss

Gen3_benchmark = DynGenerator(2, #number
                                :Gen2, #name
                                nodes_benchmark[3],
                                1.0, #ω_re
                                1.0059, #V_ref
                                1.0, #P_ref
                                Mach3_benchmark, #machine
                                Shaft3_benchmark, #shaft
                                gen3_avr_benchmark_m, #avr
                                fixed_tg, #tg
                                no_pss) #pss

Gen2_benchmark_fixedavr = DynGenerator(1, #number
                              :Gen1, #name
                              nodes_benchmark[2],
                              1.0, #ω_ref
                              1.0142, #V_ref
                              1.0, #P_ref
                              Mach2_benchmark, #machine
                              Shaft2_benchmark, #shaft
                              gen2_avr_fixed, #avr
                              fixed_tg, #tg
                              no_pss) #pss

Gen3_benchmark_fixedavr = DynGenerator(2, #number
                                :Gen2, #name
                                nodes_benchmark[3],
                                1.0, #ω_ref
                                1.00586, #V_ref
                                1.0, #P_ref
                                Mach3_benchmark, #machine
                                Shaft3_benchmark, #shaf
                                gen3_avr_fixed, #avr
                                fixed_tg, #tg
                                no_pss) #pss

Gen_bench = DynGenerator(1, #number
                              :Gen1, #name
                              nodes_bench[2],
                              1.0, #ω_ref
                              1.0, #V_ref
                              1.0, #P_ref
                              Mach2_benchmark, #machine
                              Shaft2_benchmark, #shaft
                              gen_avr_bench, #avr
                              fixed_tg, #tg
                              no_pss) #pss
                                gen3_avr_benchmark_m, #avr
                                fixed_tg, #tg
                                no_pss) #pss

################### Sources #####################

inf_gen = StaticSource(1, #number
                 :InfBus, #name
                 nodes_OMIB[1],#bus
                 1.05, #VR
                 0.0, #VI
                 0.000005) #Xth

inf_gen_kundur = StaticSource(1,
                        :InfBus,
                        nodes_kundur[2],
                        1.01,
                        0.0,
                        0.000005)

inf_gen_benchmark = StaticSource(1, #number
                 :InfBus, #name
                 nodes_benchmark[1],#bus
                 1.02, #VR
                 0.0, #VI
                 0.000001) #Xth


################### Sources #####################
inf_gen_kundur = Source(1,
                        :InfBus,
                        nodes_kundur[2],
                        1.01,
                        0.0,
                        0.000005)

inf_gen_benchmark = Source(1, #number
                 :InfBus, #name
                 nodes_benchmark[1],#bus
                 1.02, #VR
                 0.0, #VI
                 0.000001) #Xth


################### Sources #####################
OMIB =  DynamicSystem(nodes_OMIB, branch_OMIB, [Gen1AVR], [inf_gen], 100.0, 60.0)
OMIBA = DynamicSystem(nodes_OMIB, branch_OMIB, [Gen5AVR], [inf_gen], 100.0, 60.0)
OMIB_Kundur = DynamicSystem(nodes_kundur, branch_kundur, [GenKundur], [inf_gen_kundur], 100.0, 60.0)
OMIB_FullKundur = DynamicSystem(nodes_kundur, branch_kundur, [GenFullKundur], [inf_gen_kundur], 100.0, 60.0
TwoBus_Benchmark = DynamicSystem(nodes_bench, #buses
                                   branch_bench, #branches
                                   [Gen_bench], #dynamic injections
                                   [inf_gen_benchmark, loads_bench[1]], #static injections
                                   100.0, #Sbase
                                   60.0) #fbase

ThreeBus_Benchmark = DynamicSystem(nodes_benchmark, #buses
                                  branch_benchmark, #branches
                                  [Gen2_benchmark, Gen3_benchmark], #dynamic injections
                                  [inf_gen_benchmark, loads_benchmark[1], loads_benchmark[2], loads_benchmark[3]], #static injections
                                  100.0, #Sbase
                                  60.0) #fbase

ThreeBus_Benchmark = DynamicSystem(nodes_benchmark, #buses
                                  branch_benchmark, #branches
                                  [Gen2_benchmark, Gen3_benchmark], #dynamic injections
                                  [inf_gen_benchmark, loads_benchmark[1], loads_benchmark[2], loads_benchmark[3]], #static injections
                                  100.0, #Sbase
                                  60.0) #fbase

ThreeBus_Benchmark_tg = DynamicSystem(nodes_benchmark, #buses
                                  branch_benchmark, #branches
                                  [Gen2_benchmark_tg, Gen3_benchmark], #dynamic injections
                                  [inf_gen_benchmark, loads_benchmark[1], loads_benchmark[2], loads_benchmark[3]], #static injections
                                  100.0, #Sbase
                                  60.0) #fbase

ThreeBus_Benchmark_fixedavr = DynamicSystem(nodes_benchmark, #buses
                                branch_benchmark, #branches
                                [Gen2_benchmark_fixedavr, Gen3_benchmark_fixedavr], #dynamic injections
                                [inf_gen_benchmark,loads_benchmark[1], loads_benchmark[2], loads_benchmark[3]], #static injections
                                100.0, #Sbase
                                60.0) #fbase
ThreeBus_Benchmark = DynamicSystem(nodes_benchmark, #buses
                                   branch_benchmark, #branches
                                   [Gen2_benchmark, Gen3_benchmark], #dynamic injections
                                   [inf_gen_benchmark, loads_benchmark[1], loads_benchmark[2], loads_benchmark[3]], #static injections
                                   100.0, #Sbase
                                   60.0) #fbase
