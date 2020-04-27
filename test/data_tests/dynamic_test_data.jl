using PowerSystems

######################################
############ Generators ##############
######################################

######## Machine Data #########

machine_OMIB() = BaseMachine(
    0.0, #R
    0.2995, #Xd_p
    0.7087, #eq_p
    100.0,
)  #MVABase

machine_4th() = OneDOneQMachine(
    0.0, #R
    1.3125, #Xd
    1.2578, #Xq
    0.1813, #Xd_p
    0.25, #Xq_p
    5.89, #Td0_p
    0.6, #Tq0_p
    100.0,
)   #MVABase

machine_6th() = SimpleMarconatoMachine(
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

machine_8th() = MarconatoMachine(
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

machine_anderson() = AndersonFouadMachine(
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

machine_kundur() = SimpleFullMachine(
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

machine_full_kundur() = FullMachine(
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

shaft_damping() = SingleMass(
    3.148, #H
    2.0,
) #D

shaft_no_damping() = SingleMass(
    3.01, #H (M = 6.02 -> H = M/2)
    0.0,
) #D

shaft_fivemass() = FiveMassShaft(
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

pss_none() = PSSFixed(0.0)

######## TG Data #########

tg_none() = TGFixed(1.0) #eff

tg_type1() = TGTypeI(
    0.02, #R
    0.1, #Ts
    0.45, #Tc
    0.0, #T3
    0.0, #T4
    50.0, #T5
    0.3, #P_min
    1.2,
) #P_max

tg_type2() = TGTypeII(
    0.05, #R
    2.0, #T1
    1.0, #T2
    1.5, #τ_max
    0.1,
) #τ_min

########  AVR Data #########

avr_none() = AVRFixed(0.0)

avr_propr() = AVRSimple(5000.0) #Kv

avr_fixed() = AVRFixed(1.05) #Emf

avr_type1() = AVRTypeI(
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

avr_type2() = AVRTypeII(
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

######################################
############# Inverters ##############
######################################

###### Converter Data ######
converter_low_power() = AverageConverter(
    v_rated = 690.0,
    s_rated = 2.75,
)

converter_high_power() = AverageConverter(
    v_rated = 138.0,
    s_rated = 100.0,
)

###### DC Source Data ######
dc_source_lv() = FixedDCSource(voltage = 600.0) #Not in the original data, guessed.
dc_source_hv() = FixedDCSource(voltage = 1500.0) #Not in the original data, guessed.

###### Filter Data ######
filter() = LCLFilter(
    lf = 0.08,
    rf = 0.003,
    cf = 0.074,
    lg = 0.2,
    rg = 0.01,
)

###### PLL Data ######
pll() = KauraPLL(
    ω_lp = 500.0, #Cut-off frequency for LowPass filter of PLL filter.
    k_p = 0.084,  #PLL proportional gain
    k_i = 4.69,   #PLL integral gain
)

###### Outer Control ######
function outer_control()
    function virtual_inertia()
        return VirtualInertia(
        Ta = 2.0,
        kd = 400.0,
        kω = 20.0,
        ωb = 2 * pi * 50.0,
    )
    end
    function reactive_droop()
        return ReactivePowerDroop(
        kq = 0.2,
        ωf = 1000.0,
    )
    end
    return OuterControl(virtual_inertia(), reactive_droop())
end

######## Inner Control ######
inner_control() = CurrentControl(
    kpv = 0.59,     #Voltage controller proportional gain
    kiv = 736.0,    #Voltage controller integral gain
    kffv = 0.0,     #Binary variable enabling the voltage feed-forward in output of current controllers
    rv = 0.0,       #Virtual resistance in pu
    lv = 0.2,       #Virtual inductance in pu
    kpc = 1.27,     #Current controller proportional gain
    kic = 14.3,     #Current controller integral gain
    kffi = 0.0,     #Binary variable enabling the current feed-forward in output of current controllers
    ωad = 50.0,     #Active damping low pass filter cut-off frequency
    kad = 0.2,
)
