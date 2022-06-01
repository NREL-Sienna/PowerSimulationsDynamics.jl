using PowerSystems

######################################
############ Generators ##############
######################################

######## Machine Data #########

machine_classic() = BaseMachine(
    0.0, #R
    0.2995, #Xd_p
    0.7087,
)  #MVABase

machine_oneDoneQ() = OneDOneQMachine(
    0.0, #R
    1.3125, #Xd
    1.2578, #Xq
    0.1813, #Xd_p
    0.25, #Xq_p
    5.89, #Td0_p
    0.6, #Tq0_p
)

machine_simple_marconato() = SimpleMarconatoMachine(
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
    0.0,
) #MVABase

machine_marconato() = MarconatoMachine(
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
    0.0,
) #MVABase

machine_anderson() = AndersonFouadMachine(
    0.0, #R
    0.8979, #Xd
    0.646, #Xq
    0.2995, #Xd_p
    0.646, #Xq_p
    0.23, #Xd_pp
    0.4, #Xq_pp
    3.0, #Td0_p
    0.1, #Tq0_p
    0.01, #Td0_pp
    0.033, #Tq0_pp
)

machine_simple_anderson() = SimpleAFMachine(
    0.0, #R
    0.8979, #Xd
    0.646, #Xq
    0.2995, #Xd_p
    0.646, #Xq_p
    0.23, #Xd_pp
    0.4, #Xq_pp
    3.0, #Td0_p
    0.1, #Tq0_p
    0.01, #Td0_pp
    0.033, #Tq0_pp
)

machine_genrou() = RoundRotorExponential(
    R = 0.0,
    Td0_p = 8.0,
    Td0_pp = 0.03,
    Tq0_p = 0.4,
    Tq0_pp = 0.05,
    Xd = 1.8,
    Xq = 1.7,
    Xd_p = 0.3,
    Xq_p = 0.55,
    Xd_pp = 0.25,
    Xl = 0.2,
    Se = (0.0, 0.0),
)

#Not available yet
#=
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
    #0.0006, #R_f
    0.003, #R_f
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
)
=#

machine_multi_ref() = BaseMachine(
    0.0, #R
    0.2995, #Xd_p
    1.0901,
)  #MVABase

machine_multi() = BaseMachine(
    0.0, #R
    0.2995, #Xd_p
    0.9516,
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

shaft_genrou() = SingleMass(H = 6.175, D = 0.05)

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

pss_stab1() = PSY.STAB1(
    KT = 6.0,
    T = 1.5,
    T1T3 = 13.3,
    T3 = 0.0447,
    T2T4 = 13.3,
    T4 = 0.0447,
    H_lim = 0.2
)

######## TG Data #########

tg_none() = TGFixed(1.0) #eff

tg_type1() = TGTypeI(
    0.02, #R
    0.1, #Ts
    0.45, #Tc
    0.0, #T3
    12.0, #T4
    50.0, #T5
    (min = 0.3, max = 1.2), #P_lims
)

tg_type2() = TGTypeII(
    0.05, #R
    1.0, #T1
    2.0, #T2
    (min = 0.1, max = 1.5), #τ_lims
)

########  AVR Data #########

avr_none() = AVRFixed(0.0)

avr_propr() = AVRSimple(500.0) #Kv

avr_fixed() = AVRFixed(1.05) #Emf

avr_type1() = AVRTypeI(
    20.0, #Ka - Gain
    0.01, #Ke
    0.063, #Kf
    0.2, #Ta
    0.314, #Te
    0.35, #Tf
    0.001, #Tr
    (min = -5.0, max = 5.0),
    0.0039, #Ae - 1st ceiling coefficient
    1.555,
) #Be - 2nd ceiling coefficient

avr_type2() = AVRTypeII(
    200.0, #K0 - Gain
    4.0, #T1 - 1st pole
    1.0, #T2 - 1st zero
    0.006, #T3 - 2nd pole
    0.06, #T4 - 2nd zero
    0.0001, #Te - Field current time constant
    0.0001, #Tr - Measurement time constant
    (min = -50.0, max = 50.0),
    0.0, #Ae - 1st ceiling coefficient
    0.0,
) #Be - 2nd ceiling coefficient

avr_exst1() = PSY.EXST1(
    Tr = 0.01,
    Vi_lim = (-5.0, 5.0),
    Tc = 10.0,
    Tb = 20.0,
    Ka = 200.0,
    Ta = 0.1,
    Vr_lim = (0.0, 6.0),
    Kc = 0.0,
    Kf = 0.0,
    Tf = 0.1,
)

######################################
############# Inverters ##############
######################################

###### Converter Data ######
converter_low_power() = AverageConverter(rated_voltage = 690.0, rated_current = 2.75)

converter_high_power() = AverageConverter(rated_voltage = 138.0, rated_current = 100.0)

###### DC Source Data ######
dc_source_lv() = FixedDCSource(voltage = 600.0) #Not in the original data, guessed.
dc_source_hv() = FixedDCSource(voltage = 1500.0) #Not in the original data, guessed.

###### Filter Data ######
filt() = LCLFilter(lf = 0.08, rf = 0.003, cf = 0.074, lg = 0.2, rg = 0.01)
filt_gfoll() = LCLFilter(lf = 0.009, rf = 0.016, cf = 2.5, lg = 0.002, rg = 0.003)

###### PLL Data ######
pll() = KauraPLL(
    ω_lp = 500.0, #Cut-off frequency for LowPass filter of PLL filter.
    kp_pll = 0.084,  #PLL proportional gain
    ki_pll = 4.69,   #PLL integral gain
)

reduced_pll() = ReducedOrderPLL(
    ω_lp = 1.32 * 2 * pi * 50, #Cut-off frequency for LowPass filter of PLL filter.
    kp_pll = 2.0,  #PLL proportional gain
    ki_pll = 20.0,   #PLL integral gain
)

no_pll() = PSY.FixedFrequency()

###### Outer Control ######
function outer_control()
    function virtual_inertia()
        return VirtualInertia(Ta = 2.0, kd = 400.0, kω = 20.0)
    end
    function reactive_droop()
        return ReactivePowerDroop(kq = 0.2, ωf = 1000.0)
    end
    return OuterControl(virtual_inertia(), reactive_droop())
end

function outer_control_nopll()
    function virtual_inertia()
        return VirtualInertia(Ta = 2.0, kd = 0.0, kω = 20.0)
    end
    function reactive_droop()
        return ReactivePowerDroop(kq = 0.2, ωf = 1000.0)
    end
    return OuterControl(virtual_inertia(), reactive_droop())
end

function outer_control_droop()
    function active_droop()
        return PSY.ActivePowerDroop(Rp = 0.05, ωz = 2 * pi * 5)
    end
    function reactive_droop()
        return ReactivePowerDroop(kq = 0.2, ωf = 1000.0)
    end
    return OuterControl(active_droop(), reactive_droop())
end

function outer_control_gfoll()
    function active_pi()
        return ActivePowerPI(Kp_p = 2.0, Ki_p = 30.0, ωz = 0.132 * 2 * pi * 50)
    end
    function reactive_pi()
        return ReactivePowerPI(Kp_q = 2.0, Ki_q = 30.0, ωf = 0.132 * 2 * pi * 50.0)
    end
    return OuterControl(active_pi(), reactive_pi())
end

######## Inner Control ######
inner_control() = VoltageModeControl(
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

current_mode_inner() = CurrentModeControl(
    kpc = 0.37,     #Current controller proportional gain
    kic = 0.7,     #Current controller integral gain
    kffv = 1.0,     #Binary variable enabling the voltage feed-forward in output of current controllers
)

#########################################
####### Generic Renewable Models ########
#########################################

function outer_control_TypeAB()
    function active_ab()
        return ActiveRenewableControllerAB(
            bus_control = 0,
            from_branch_control = 0,
            to_branch_control = 0,
            branch_id_control = "0",
            Freq_Flag = 0,
            K_pg = 1.0,
            K_ig = 0.05,
            T_p = 0.25,
            fdbd_pnts = (-0.01, 0.01),
            fe_lim = (-99.0, 99.0),
            P_lim = (0.0, 1.2),
            T_g = 0.1,
            D_dn = 0.0,
            D_up = 0.0,
            dP_lim = (-99.0, 99.0),
            P_lim_inner = (0.0, 1.2),
            T_pord = 0.02,
        )
    end
    function reactive_ab()
        return ReactiveRenewableControllerAB(
            bus_control = 0,
            from_branch_control = 0,
            to_branch_control = 0,
            branch_id_control = "0",
            VC_Flag = 0,
            Ref_Flag = 0,
            PF_Flag = 0,
            V_Flag = 0,
            T_fltr = 0.02,
            K_p = 18.0,
            K_i = 50.0,
            T_ft = 0.0,
            T_fv = 0.05,
            V_frz = 0.0,
            R_c = 0.0,
            X_c = 0.0,
            K_c = 0.0,
            e_lim = (-0.1, 0.1),
            dbd_pnts = (-0.01, 0.01),
            Q_lim = (-0.43, 0.43),
            T_p = 0.0,
            Q_lim_inner = (-0.44, 0.44),
            V_lim = (0.9, 1.05),
            K_qp = 0.0,
            K_qi = 0.01,
        )
    end
    return OuterControl(active_ab(), reactive_ab())
end

function outer_control_TypeAB_FreqFlag()
    function active_ab_fflag()
        return ActiveRenewableControllerAB(
            bus_control = 0,
            from_branch_control = 0,
            to_branch_control = 0,
            branch_id_control = "0",
            Freq_Flag = 1,
            K_pg = 1.0,
            K_ig = 0.05,
            T_p = 0.25,
            fdbd_pnts = (-0.01, 0.01),
            fe_lim = (-99.0, 99.0),
            P_lim = (0.0, 1.2),
            T_g = 0.1,
            D_dn = 0.0,
            D_up = 0.0,
            dP_lim = (-99.0, 99.0),
            P_lim_inner = (0.0, 1.2),
            T_pord = 0.02,
        )
    end
    function reactive_ab()
        return ReactiveRenewableControllerAB(
            bus_control = 0,
            from_branch_control = 0,
            to_branch_control = 0,
            branch_id_control = "0",
            VC_Flag = 0,
            Ref_Flag = 0,
            PF_Flag = 0,
            V_Flag = 0,
            T_fltr = 0.02,
            K_p = 18.0,
            K_i = 50.0,
            T_ft = 0.0,
            T_fv = 0.05,
            V_frz = 0.0,
            R_c = 0.0,
            X_c = 0.0,
            K_c = 0.0,
            e_lim = (-0.1, 0.1),
            dbd_pnts = (-0.01, 0.01),
            Q_lim = (-0.43, 0.43),
            T_p = 0.0,
            Q_lim_inner = (-0.44, 0.44),
            V_lim = (0.9, 1.05),
            K_qp = 0.0,
            K_qi = 0.01,
        )
    end
    return OuterControl(active_ab_fflag(), reactive_ab())
end

inner_ctrl_typeB() = RECurrentControlB(
    Q_Flag = 0,
    PQ_Flag = 0,
    Vdip_lim = (-99.0, 99.0),
    T_rv = 0.02,
    dbd_pnts = (-0.01, 0.01),
    K_qv = 1.0,
    Iqinj_lim = (-1.1, 1.1),
    V_ref0 = 0.0,
    K_vp = 10.0,
    K_vi = 60.0,
    T_iq = 0.02,
    I_max = 1.11,
)

converter_regca() = RenewableEnergyConverterTypeA(
    T_g = 0.02,
    Rrpwr = 10.0,
    Brkpt = 0.9,
    Zerox = 0.4,
    Lvpl1 = 1.22,
    Vo_lim = 1.2,
    Lv_pnts = (0.5, 0.9),
    Io_lim = -1.3,
    T_fltr = 0.2,
    K_hv = 0.0,
    Iqr_lims = (-100.0, 100.0),
    Accel = 0.7,
    Lvpl_sw = 0,
    R_source = 0.000,
    X_source = 1.0e5,
)

filt_current() = RLFilter(rf = 0.0, lf = 0.0)

####### Loads ########

Ind_Motor(load) = PSY.SingleCageInductionMachine(
    name = PSY.get_name(load),
    R_s = 0.013,
    R_r = 0.009,
    X_ls = 0.067,
    X_lr = 0.17,
    X_m = 3.8,
    H = 1.5,
    A = 1.0,
    B = 0.0,
    base_power = 1000.0,
)

####### Devices #######

function dyn_gen_second_order(generator)
    return PSY.DynamicGenerator(
        name = get_name(generator),
        ω_ref = 1.0, # ω_ref,
        machine = machine_oneDoneQ(), #machine
        shaft = shaft_no_damping(), #shaft
        avr = avr_type1(), #avr
        prime_mover = tg_none(), #tg
        pss = pss_none(), #pss
    )
end

function inv_case78(static_device)
    return DynamicInverter(
        name = get_name(static_device),
        ω_ref = 1.0, # ω_ref,
        converter = converter_high_power(), #converter
        outer_control = outer_control(), #outer control
        inner_control = inner_control(), #inner control voltage source
        dc_source = dc_source_lv(), #dc source
        freq_estimator = pll(), #pll
        filter = filt(), #filter
    )
end
