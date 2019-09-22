mach_test = BaseMachine(   0.0, #R
                       0.2995, #Xd_p
                       0.588, #eq_p
                       100.0)  #MVABase

mach2_test = OneDOneQMachine(0.0, #R
                      0.8, #Xd
                       0.2, #Xq
                       0.1813, #Xd_p
                       0.05, #Xq_p
                       5.89, #Td0_p
                       0.6, #Tq0_p
                       100.0)   #MVABase

shaft_test = SingleMass(3.148, #H
                       2.0) #D

shaft_test2 = SingleMass(3.01, #H
                      2.0) #D

pss_test = PSSFixed(0.0)
tg_test = TGFixed(1.0)
avr_test = AVRFixed(1.0)
avr_test2 = AVRFixed(0.26854)
#avr_test3 = AVRFixed(0.3337)
avr_test3 = AVRFixed(0.34992)

load_test = PowerLoad("Bus3", true, nodes_test[1], PowerSystems.ConstantPower, 0.0, 0.0, 0.0, 0.0)
load_test2 = PowerLoad("Bus3", true, nodes_test2[3], PowerSystems.ConstantPower, 0.6, 0.05, 0.6, 0.05)

load_test3a = PowerLoad("Bus3", true, nodes_test2[2], PowerSystems.ConstantPower, 0.3, 0.05, 0.6, 0.05)
load_test3b = PowerLoad("Bus3", true, nodes_test2[3], PowerSystems.ConstantPower, 0.3, 0.05, 0.6, 0.05)

Gen_test = DynGenerator(1, #Number
                 :TestGen,
                 nodes_test[2], #bus
                 1.0, # ω_ref,
                 1.0, #V_ref
                 1.0, #P_ref
                 mach_test,
                 shaft_test,
                 avr_test, #avr
                 tg_test, #tg
                 pss_test)

Gen_test2 = DynGenerator(1, #Number
                :TestGen,
                nodes_test[2], #bus
                1.0, # ω_ref,
                1.0, #V_ref
                1.0, #P_ref
                mach2_test,
                shaft_test,
                avr_test2, #avr
                tg_test, #tg
                pss_test)

Gen_test3 = DynGenerator(1, #Number
              :TestGen,
              nodes_test[2], #bus
              1.0, # ω_ref,
              1.0, #V_ref
              0.5, #P_ref
              mach2_test,
              shaft_test,
              avr_test3, #avr
              tg_test, #tg
              pss_test)

OMIB_test =  DynamicSystem(nodes_test, branch_test, [Gen_test], [load_test], 100.0, 60.0)
OMIB_test2 =  DynamicSystem(nodes_test, branch_test, [Gen_test2], [load_test], 100.0, 60.0)
OMIB_test_load = DynamicSystem(nodes_test2, branch_test2, [Gen_test3], [load_test2], 100.0, 60.0)
OMIB_test_load2 = DynamicSystem(nodes_test2, branch_test2, [Gen_test3], [load_test3a, load_test3b], 100.0, 60.0)
