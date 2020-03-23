using PowerSystems
const PSY = PowerSystems

######################################
############ Generators ##############
######################################

######## Machine Data #########

machine_OMIB() = PSY.BaseMachine(
    0.0, #R
    0.2995, #Xd_p
    0.7087, #eq_p
    100.0,
)  #MVABase

machine_4th() = PSY.OneDOneQMachine(
    0.0, #R
    1.3125, #Xd
    1.2578, #Xq
    0.1813, #Xd_p
    0.25, #Xq_p
    5.89, #Td0_p
    0.6, #Tq0_p
    100.0,
)   #MVABase

machine_6th() = PSY.SimpleMarconatoMachine(
    0.0,
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
    100.0,
) #MVABase

machine_8th() = PSY.MarconatoMachine(
    0.0,
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
    100.0,
) #MVABase

machine_anderson() = PSY.AndersonFouadMachine(
    0.0, #R
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
    615.0,
) #MVABase

machine_kundur() = PSY.SimpleFullMachine(
    0.003, #R on Example 3.1 and 4.1 of Kundur
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
    555.0,
) #MVABase

machine_full_kundur() = PSY.FullMachine(
    0.003, #R on Example 3.1 and 4.1 of Kundur
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
    555.0,
) #MVABase

machine_multi_ref() = BaseMachine(
    0.0, #R
    0.2995, #Xd_p
    1.0901, #eq_p
    100.0,
)  #MVABase

machine_multi() = BaseMachine(
    0.0, #R
    0.2995, #Xd_p
    0.9516, #eq_p
    100.0,
)  #MVABase

######## Shaft Data #########

shaft_damping() = PSY.SingleMass(
    3.148, #H
    2.0,
) #D

shaft_no_damping() = PSY.SingleMass(
    3.01, #H (M = 6.02 -> H = M/2)
    0.0,
) #D

shaft_fivemass() = PSY.FiveMassShaft(
    3.01, #5.148, #H
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
    21.984,
) #K_ex

######## PSS Data #########

pss_none() = PSY.PSSFixed(0.0)

######## TG Data #########

tg_none() = PSY.TGFixed(1.0) #eff

tg_type1() = PSY.TGTypeI(
    0.02, #R
    0.1, #Ts
    0.45, #Tc
    0.0, #T3
    0.0, #T4
    50.0, #T5
    0.3, #P_min
    1.2,
) #P_max

tg_type2() = PSY.TGTypeII(
    0.05, #R
    2.0, #T1
    1.0, #T2
    1.5, #τ_max
    0.1,
) #τ_min

########  AVR Data #########

avr_none() = PSY.AVRFixed(0.0)

avr_propr() = PSY.AVRSimple(5000.0) #Kv

avr_fixed() = PSY.AVRFixed(1.05) #Emf

avr_type1() = PSY.AVRTypeI(
    20.0, #Ka - Gain
    0.01, #Ke
    0.063, #Kf
    0.2, #Ta
    0.314, #Te
    0.35, #Tf
    0.001, #Tr
    5.0, #Vrmax
    -5.0, #Vrmin
    0.0039, #Ae - 1st ceiling coefficient
    1.555,
) #Be - 2nd ceiling coefficient

avr_type2() = PSY.AVRTypeII(
    20.0, #K0 - Gain
    0.2, #T1 - 1st pole
    0.063, #T2 - 1st zero
    0.35, #T3 - 2nd pole
    0.01, #T4 - 2nd zero
    0.314, #Te - Field current time constant
    0.001, #Tr - Measurement time constant
    5.0, #Vrmax
    -5.0, #Vrmin
    0.0039, #Ae - 1st ceiling coefficient
    1.555,
) #Be - 2nd ceiling coefficient

###### Dynamic Generators constructors ######

function dyn_gen_OMIB(nodes)
    return PSY.DynamicGenerator(
        1, #Number
        "Case1Gen",
        nodes[2], #bus
        1.0, # ω_ref,
        1.0, #V_ref
        0.5, #P_ref
        0.0, #Q_ref
        machine_OMIB(), #machine
        shaft_damping(), #shaft
        avr_none(), #avr
        tg_none(), #tg
        pss_none(),
    ) #pss
end

function dyn_gen2_case2(nodes)
    return PSY.DynamicGenerator(
        1, #Number
        "Case2Gen2",
        nodes[2], #bus
        1.0, # ω_ref,
        1.0142, #V_ref
        1.0, #P_ref
        0.0, #Q_ref
        machine_4th(), #machine
        shaft_no_damping(), #shaft
        avr_type1(), #avr
        tg_none(), #tg
        pss_none(),
    ) #pss
end

function dyn_gen3_case2(nodes)
    return PSY.DynamicGenerator(
        2, #Number
        "Case2Gen3",
        nodes[3], #bus
        1.0, # ω_ref,
        1.0059, #V_ref
        1.0, #P_ref
        0.0, #Q_ref
        machine_4th(), #machine
        shaft_no_damping(), #shaft
        avr_type1(), #avr
        tg_none(), #tg
        pss_none(),
    ) #pss
end

function dyn_gen2_case3(nodes)
    return PSY.DynamicGenerator(
        1, #Number
        "Case3Gen2",
        nodes[2], #bus
        1.0, # ω_ref,
        1.0142, #V_ref
        1.0, #P_ref
        0.0, #Q_ref
        machine_6th(), #machine
        shaft_no_damping(), #shaft
        avr_type1(), #avr
        tg_none(), #tg
        pss_none(),
    ) #pss
end

function dyn_gen3_case3(nodes)
    return PSY.DynamicGenerator(
        2, #Number
        "Case3Gen3",
        nodes[3], #bus
        1.0, # ω_ref,
        1.0059, #V_ref
        1.0, #P_ref
        0.0, #Q_ref
        machine_6th(), #machine
        shaft_no_damping(), #shaft
        avr_type1(), #avr
        tg_none(), #tg
        pss_none(),
    ) #pss
end

function dyn_gen2_case4(nodes)
    return PSY.DynamicGenerator(
        1, #Number
        "Case4Gen2",
        nodes[2], #bus
        1.0, # ω_ref,
        1.0142, #V_ref
        1.0, #P_ref
        0.0, #Q_ref
        machine_8th(), #machine
        shaft_no_damping(), #shaft
        avr_type1(), #avr
        tg_none(), #tg
        pss_none(),
    ) #pss
end

function dyn_gen3_case4(nodes)
    return PSY.DynamicGenerator(
        2, #Number
        "Case4Gen3",
        nodes[3], #bus
        1.0, # ω_ref,
        1.0059, #V_ref
        1.0, #P_ref
        0.0, #Q_ref
        machine_8th(), #machine
        shaft_no_damping(), #shaft
        avr_type1(), #avr
        tg_none(), #tg
        pss_none(),
    ) #pss
end

function dyn_gen_case5(nodes)
    return PSY.DynamicGenerator(
        1, #Number
        "Case5Gen",
        nodes[2], #bus
        1.0, # ω_ref,
        1.0155, #V_ref
        0.5, #P_ref
        0.0, #Q_ref
        machine_4th(), #machine
        shaft_fivemass(), #shaft
        avr_type1(), #avr
        tg_type2(), #tg
        pss_none(),
    ) #pss
end

function dyn_gen_case7(nodes)
    return PSY.DynamicGenerator(
        1, #Number
        "Case7Gen",
        nodes[2], #bus
        1.0, # ω_ref,
        1.0142, #V_ref
        0.6, #P_ref
        0.0, #Q_ref
        machine_4th(), #machine
        shaft_no_damping(), #shaft
        avr_type1(), #avr
        tg_none(), #tg
        pss_none(),
    ) #pss
end

function dyn_gen_case8(nodes)
    return PSY.DynamicGenerator(
        1, #Number
        "Case8Gen",
        nodes[2], #bus
        1.0, # ω_ref,
        1.0124, #V_ref
        0.6, #P_ref
        0.0, #Q_ref
        machine_4th(), #machine
        shaft_no_damping(), #shaft
        avr_type1(), #avr
        tg_none(), #tg
        pss_none(),
    ) #pss
end

function dyn_gen_case9(nodes)
    return PSY.DynamicGenerator(
        1, #Number
        "Case9Gen",
        nodes[2], #bus
        1.0, # ω_ref,
        1.0124, #V_ref
        0.6, #P_ref
        0.0, #Q_ref
        machine_4th(), #machine
        shaft_no_damping(), #shaft
        avr_type1(), #avr
        tg_none(), #tg
        pss_none(),
    ) #pss
end

function dyn_gen_case10_ref(nodes)
    return PSY.DynamicGenerator(
        1, #Number
        "Case10Gen1",
        nodes[1], #bus
        1.0, # ω_ref,
        1.0, #V_ref
        0.41, #P_ref
        0.0, #Q_ref
        machine_multi_ref(), #machine
        shaft_damping(), #shaft
        avr_none(), #avr
        tg_none(), #tg
        pss_none(),
    ) #pss
end

function dyn_gen_case10_2(nodes)
    return PSY.DynamicGenerator(
        1, #Number
        "Case10Gen2",
        nodes[3], #bus
        1.0, # ω_ref,
        1.0, #V_ref
        0.4, #P_ref
        0.0, #Q_ref
        machine_multi(), #machine
        shaft_damping(), #shaft
        avr_none(), #avr
        tg_type2(), #tg
        pss_none(),
    ) #pss
end

######################################
############# Inverters ##############
######################################

###### Converter Data ######

converter_DAIB() = PSY.AvgCnvFixedDC(
    690.0, #Rated Voltage
    2.75,
) #Rated MVA

converter_case78() = PSY.AvgCnvFixedDC(
    138.0, #Rated Voltage
    100.0,
) #Rated MVA

###### DC Source Data ######

dc_source_DAIB() = PSY.FixedDCSource(600.0) #Not in the original data, guessed.

dc_source_case78() = PSY.FixedDCSource(1500.0) #Not in the original data, guessed.

###### Filter Data ######

filter_test() = PSY.LCLFilter(
    0.08, #Series inductance lf in pu
    0.003, #Series resitance rf in pu
    0.074, #Shunt capacitance cf in pu
    0.2, #Series ractance rg to grid connection (#Step up transformer or similar)
    0.01,
) #Series resistance lg to grid connection (#Step up transformer or similar)

###### PLL Data ######

pll_test() = PSY.PLL(
    500.0, #ω_lp: Cut-off frequency for LowPass filter of PLL filter.
    0.084, #k_p: PLL proportional gain
    4.69,
) #k_i: PLL integral gain

###### Outer Control ######

virtual_inertia_test() = PSY.VirtualInertia(
    2.0, #Ta:: VSM inertia constant
    400.0, #kd:: VSM damping coefficient
    20.0, #kω:: Frequency droop gain in pu
    2 * pi * 50.0,
) #ωb:: Rated angular frequency

reactive_droop_test() = PSY.ReactivePowerDroop(
    0.2, #kq:: Reactive power droop gain in pu
    1000.0,
) #ωf:: Reactive power cut-off low pass filter frequency

outer_control_test() =
    PSY.VirtualInertiaQdroop(virtual_inertia_test(), reactive_droop_test())

######## Inner Control ######

vsc_test() = PSY.CombinedVIwithVZ(
    0.59, #kpv:: Voltage controller proportional gain
    736.0, #kiv:: Voltage controller integral gain
    0.0, #kffv:: Binary variable enabling the voltage feed-forward in output of current controllers
    0.0, #rv:: Virtual resistance in pu
    0.2, #lv: Virtual inductance in pu
    1.27, #kpc:: Current controller proportional gain
    14.3, #kiv:: Current controller integral gain
    0.0, #kffi:: Binary variable enabling the current feed-forward in output of current controllers
    50.0, #ωad:: Active damping low pass filter cut-off frequency
    0.2,
) #kad:: Active damping gain

###### Inverters constructors ######

function inv_DAIB(nodes)
    return PSY.DynamicInverter(
        1, #Number
        "DARCO", #name
        nodes[2], #bus
        1.0, # ω_ref,
        1.02, #V_ref
        0.5, #P_ref
        0.0, #Q_ref
        2.75, #MVABase
        converter_DAIB(), #converter
        outer_control_test(), #outer control
        vsc_test(), #inner control voltage source
        dc_source_DAIB(), #dc source
        pll_test(), #pll
        filter_test(),
    ) #filter
end

function inv_case78(nodes)
    return PSY.DynamicInverter(
        1, #Number
        "DARCO", #name
        nodes[2], #bus
        1.0, # ω_ref,
        1.02, #V_ref
        0.5, #P_ref
        0.0, #Q_ref
        2.75, #MVABase
        converter_case78(), #converter
        outer_control_test(), #outer control
        vsc_test(), #inner control voltage source
        dc_source_case78(), #dc source
        pll_test(), #pll
        filter_test(),
    ) #filter
end

function inv_case9(nodes)
    return PSY.DynamicInverter(
        1, #Number
        "DARCO", #name
        nodes[3], #bus
        1.0, # ω_ref,
        0.8, #V_ref
        0.5, #P_ref
        -0.3, #Q_ref
        100.0, #MVABase
        converter_case78(), #converter
        outer_control_test(), #outer control
        vsc_test(), #inner control voltage source
        dc_source_case78(), #dc source
        pll_test(), #pll
        filter_test(),
    ) #filter
end

######################################
######## System Constructors #########
######################################

function system_no_inv(nodes, branches, loads, sources, gens)
    #Create system with BasePower = 100 MVA and nominal frequency 60 Hz.
    sys = PSY.System(100.0, frequency = 60.0)

    #Add buses
    for bus in nodes
        PSY.add_component!(sys, bus)
    end

    #Add lines
    for lines in branches
        PSY.add_component!(sys, lines)
    end

    #Add loads
    for load in loads
        PSY.add_component!(sys, load)
    end

    #Add infinite source
    for source in sources
        PSY.add_component!(sys, source)
    end

    #Add generator
    for gen in gens
        PSY.add_component!(sys, gen)
    end
    return sys
end

function system_no_inv_no_sources(nodes, branches, loads, gens)
    #Create system with BasePower = 100 MVA and nominal frequency 60 Hz.
    sys = PSY.System(100.0, frequency = 60.0)

    #Add buses
    for bus in nodes
        PSY.add_component!(sys, bus)
    end

    #Add lines
    for lines in branches
        PSY.add_component!(sys, lines)
    end

    #Add loads
    for load in loads
        PSY.add_component!(sys, load)
    end

    #Add generator
    for gen in gens
        PSY.add_component!(sys, gen)
    end
    return sys
end

function system_DAIB(nodes, branches, sources, invs)
    #Create system with BasePower = 100 MVA and nominal frequency 50 Hz.
    sys = PSY.System(100.0, frequency = 50.0)

    #Add buses
    for bus in nodes
        PSY.add_component!(sys, bus)
    end

    #Add lines
    for lines in branches
        PSY.add_component!(sys, lines)
    end

    #Add infinite source
    for source in sources
        PSY.add_component!(sys, source)
    end

    #Add inverters
    for inv in invs
        PSY.add_component!(sys, inv)
    end
    return sys
end

function system_50Hz(nodes, branches, loads, sources, invs, gens)
    #Create system with BasePower = 100 MVA and nominal frequency 50 Hz.
    sys = PSY.System(100.0, frequency = 50.0)

    #Add buses
    for bus in nodes
        PSY.add_component!(sys, bus)
    end

    #Add lines
    for lines in branches
        PSY.add_component!(sys, lines)
    end

    #Add loads
    for load in loads
        PSY.add_component!(sys, load)
    end

    #Add infinite source
    for source in sources
        PSY.add_component!(sys, source)
    end

    #Add inverters
    for inv in invs
        PSY.add_component!(sys, inv)
    end

    #Add generators
    for gen in gens
        PSY.add_component!(sys, gen)
    end

    return sys
end

function system_case9(nodes, branches, loads, sources, invs, gens)
    #Create system with BasePower = 100 MVA and nominal frequency 50 Hz.
    sys = PSY.System(100.0, frequency = 50.0)

    #Add buses
    for bus in nodes
        PSY.add_component!(sys, bus)
    end

    #Add lines
    for lines in branches
        PSY.add_component!(sys, lines)
    end

    #Make line 3 the dynamic line
    make_dynamic_branch!(branches[3], sys)

    #Add loads
    for load in loads
        PSY.add_component!(sys, load)
    end

    #Add infinite source
    for source in sources
        PSY.add_component!(sys, source)
    end

    #Add inverters
    for inv in invs
        PSY.add_component!(sys, inv)
    end

    #Add generators
    for gen in gens
        PSY.add_component!(sys, gen)
    end

    return sys
end
