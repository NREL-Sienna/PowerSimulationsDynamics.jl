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
                        
Shaft2_benchmark = SingleMass(3.01, #H (M = 6.02 -> H = M/2)
                       0.0) #D

Shaft3_benchmark = SingleMass(3.01, #H (M = 6.02 -> H = M/2)
                      0.0) #D
                      
                      
no_pss = PSSFixed(0.0)

fixed_tg = TGFixed(1.0) #eff

gen2_avr_benchmark_m = AVRTypeI(20.0, #Ka - Gain
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
                              
Gen3_benchmark = DynGenerator(2, #number
                                :Gen2, #name
                                nodes_benchmark[3],
                                1.0, #ω_ref
                                1.0059, #V_ref
                                1.0, #P_ref
                                Mach3_benchmark, #machine
                                Shaft3_benchmark, #shaft
                                gen3_avr_benchmark_m, #avr
                                fixed_tg, #tg
                                no_pss) #pss
                                
inf_gen_benchmark = StaticSource(1, #number
                 :InfBus, #name
                 nodes_benchmark[1], #bus
                 1.02, #VR
                 0.0, #VI
                 0.000001) #Xth
                 
                 
ThreeBus_Benchmark = DynamicSystem(nodes_benchmark, #buses
                                  branch_benchmark, #branches
                                  [Gen2_benchmark, Gen3_benchmark], #dynamic injections
                                  [inf_gen_benchmark, loads_benchmark[1], loads_benchmark[2], loads_benchmark[3]], 
                                  100.0, #Sbase
                                  60.0) #fbase
                                  