
get_params(
    d::DynamicWrapper{T},
) where {T <: Union{PSY.DynamicGenerator, PSY.DynamicInverter}} =
    vcat(get_params(get_static_device(d)), get_params(get_dynamic_device(d)))
get_params_symbol(
    d::DynamicWrapper{T},
) where {T <: Union{PSY.DynamicGenerator, PSY.DynamicInverter}} =
    vcat(:Q_ref, get_params_symbol(get_dynamic_device(d)))
get_params(d::DynamicWrapper) = get_params(get_dynamic_device(d))
get_params(d::StaticWrapper) = get_params(get_device(d))
get_n_params(x::BranchWrapper) = get_n_params(get_branch(x))
get_params(x::BranchWrapper) = get_params(PSY.get_branch(get_branch(x)))
get_n_params(x::PSY.DynamicBranch) = get_n_params(PSY.get_branch(x))
get_n_params(x::PSY.Line) = 2
get_params(x::PSY.Line) = [PSY.get_r(x), PSY.get_x(x)]

function get_params(d::StaticLoadWrapper)
    loads = get_loads(d)
    bus = PSY.get_bus(d)
    V_ref = PSY.get_magnitude(bus)
    Θ_ref = PSY.get_angle(bus)
    sys_base_power = get_system_base_power(d)
    P_power = 0.0
    P_current = 0.0
    P_impedance = 0.0
    Q_power = 0.0
    Q_current = 0.0
    Q_impedance = 0.0
    for ld in loads
        base_power_conversion = PSY.get_base_power(ld) / sys_base_power
        if isa(ld, PSY.PowerLoad)
            P_power += PSY.get_active_power(ld) * base_power_conversion
            Q_power += PSY.get_reactive_power(ld) * base_power_conversion
        elseif isa(ld, PSY.StandardLoad)
            P_impedance += PSY.get_impedance_active_power(ld) * base_power_conversion
            Q_impedance += PSY.get_impedance_reactive_power(ld) * base_power_conversion
            P_current += PSY.get_current_active_power(ld) * base_power_conversion
            Q_current += PSY.get_current_reactive_power(ld) * base_power_conversion
            P_power += PSY.get_constant_active_power(ld) * base_power_conversion
            Q_power += PSY.get_constant_reactive_power(ld) * base_power_conversion
        end
    end

    return [V_ref, Θ_ref, P_power, P_current, P_impedance, Q_power, Q_current, Q_impedance]
end

# TODO - temporary for Dynamic components that have not yet been modified to use parameters.
# Allows the 
get_n_params(x::PSY.DynamicComponent) = 0
get_params(::PSY.DynamicComponent) = Float64[]
get_params_symbol(::PSY.DynamicComponent) = Symbol[]
get_n_params(x::PSY.ActivePowerControl) = 0
get_params(::PSY.ActivePowerControl) = Float64[]
get_params_symbol(::PSY.ActivePowerControl) = Symbol[]
get_n_params(x::PSY.ReactivePowerControl) = 0
get_params(::PSY.ReactivePowerControl) = Float64[]
get_params_symbol(::PSY.ReactivePowerControl) = Symbol[]

get_n_params(g::PSY.DynamicInverter) =
    3 + get_n_params(PSY.get_converter(g)) +
    get_n_params(PSY.get_outer_control(g)) + get_n_params(PSY.get_inner_control(g)) +
    get_n_params(PSY.get_dc_source(g)) + get_n_params(PSY.get_freq_estimator(g)) +
    get_n_params(PSY.get_filter(g))

#INVERTERS 
function get_params(g::PSY.DynamicInverter)
    refs = [
        PSY.get_V_ref(g),
        PSY.get_ω_ref(g),
        PSY.get_P_ref(g),
    ]
    vcat(
        refs,
        get_params(PSY.get_converter(g)),
        get_params(PSY.get_outer_control(g)),
        get_params(PSY.get_inner_control(g)),
        get_params(PSY.get_dc_source(g)),
        get_params(PSY.get_freq_estimator(g)),
        get_params(PSY.get_filter(g)),
    )
end
get_params_symbol(g::PSY.DynamicInverter) = vcat(
    [:V_ref, :ω_ref, :P_ref],
    get_params_symbol(PSY.get_converter(g)),
    get_params_symbol(PSY.get_outer_control(g)),
    get_params_symbol(PSY.get_inner_control(g)),
    get_params_symbol(PSY.get_dc_source(g)),
    get_params_symbol(PSY.get_freq_estimator(g)),
    get_params_symbol(PSY.get_filter(g)),
)

#FILTERS 
get_n_params(x::PSY.LCLFilter) = 5
get_params(x::PSY.LCLFilter) =
    [PSY.get_lf(x), PSY.get_rf(x), PSY.get_cf(x), PSY.get_lg(x), PSY.get_rg(x)]
get_params_symbol(::PSY.LCLFilter) =
    [:lf_Filter, :rf_Filter, :cf_Filter, :lg_Filter, :rg_Filter]

get_n_params(x::PSY.RLFilter) = 2
get_params(x::PSY.RLFilter) = [PSY.get_rf(x), PSY.get_lf(x)]
get_params_symbol(::PSY.RLFilter) = [:rf_Filter, :lf_Filter]
#OUTER CONTROL
get_n_params(x::PSY.OuterControl) =
    get_n_params(PSY.get_active_power_control(x)) +
    get_n_params(PSY.get_reactive_power_control(x))
get_params(x::PSY.OuterControl) = vcat(
    get_params(PSY.get_active_power_control(x)),
    get_params(PSY.get_reactive_power_control(x)),
)
get_params_symbol(x::PSY.OuterControl) = vcat(
    get_params_symbol(PSY.get_active_power_control(x)),
    get_params_symbol(PSY.get_reactive_power_control(x)),
)
#ACTIVE POWER CONTROL
get_n_params(::PSY.VirtualInertia) = 3
get_params(x::PSY.VirtualInertia) = [PSY.get_Ta(x), PSY.get_kd(x), PSY.get_kω(x)]
get_params_symbol(::PSY.VirtualInertia) =
    [:Ta_OuterControl, :kd_OuterControl, :kω_OuterControl]

get_n_params(::PSY.ActiveRenewableControllerAB) = 17
get_params(x::PSY.ActiveRenewableControllerAB) =
    [PSY.get_K_pg(x),
        PSY.get_K_ig(x),
        PSY.get_T_p(x),
        PSY.get_fdbd_pnts(x)[1],
        PSY.get_fdbd_pnts(x)[2],
        PSY.get_fe_lim(x)[1],
        PSY.get_fe_lim(x)[2],
        PSY.get_P_lim(x)[1],
        PSY.get_P_lim(x)[2],
        PSY.get_T_g(x),
        PSY.get_D_dn(x),
        PSY.get_D_up(x),
        PSY.get_dP_lim(x)[1],
        PSY.get_dP_lim(x)[2],
        PSY.get_P_lim_inner(x)[1],
        PSY.get_P_lim_inner(x)[2],
        PSY.get_T_pord(x)]
get_params_symbol(::PSY.ActiveRenewableControllerAB) = [:K_pg_OuterControl,
    :K_ig_OuterControl,
    :T_p_ap_OuterControl, #modified to make unique
    :fdbd1_OuterControl,
    :fdbd2_OuterControl,
    :fe_min_OuterControl,
    :fe_max_OuterControl,
    :P_min_OuterControl,
    :P_max_OuterControl,
    :T_g_ap_OuterControl, #modified to make unique
    :D_dn_OuterControl,
    :D_up_OuterControl,
    :dP_min_OuterControl,
    :dP_max_OuterControl,
    :P_min_inner_OuterControl,
    :P_max_inner_OuterControl,
    :T_pord_OuterControl]

get_n_params(::PSY.ReactiveRenewableControllerAB) = 22
get_params(x::PSY.ReactiveRenewableControllerAB) = [PSY.get_T_fltr(x),
    PSY.get_K_p(x),
    PSY.get_K_i(x),
    PSY.get_T_ft(x),
    PSY.get_T_fv(x),
    PSY.get_V_frz(x),
    PSY.get_R_c(x),
    PSY.get_X_c(x),
    PSY.get_K_c(x),
    PSY.get_e_lim(x)[1],
    PSY.get_e_lim(x)[2],
    PSY.get_dbd_pnts(x)[1],
    PSY.get_dbd_pnts(x)[2],
    PSY.get_Q_lim(x)[1],
    PSY.get_Q_lim(x)[2],
    PSY.get_T_p(x),
    PSY.get_Q_lim_inner(x)[1],
    PSY.get_Q_lim_inner(x)[2],
    PSY.get_V_lim(x)[1],
    PSY.get_V_lim(x)[2],
    PSY.get_K_qp(x),
    PSY.get_K_qi(x)]
get_params_symbol(::PSY.ReactiveRenewableControllerAB) = [:T_fltr_OuterControl,
    :K_p_OuterControl,
    :K_i_OuterControl,
    :T_ft_OuterControl,
    :T_fv_OuterControl,
    :V_frz_OuterControl,
    :R_c_OuterControl,
    :X_c_OuterControl,
    :K_c_OuterControl,
    :e_min_OuterControl,
    :e_max_OuterControl,
    :dbd_pnts1_OuterControl,
    :dbd_pnts2_OuterControl,
    :Q_min_OuterControl,
    :Q_max_OuterControl,
    :T_p_OuterControl,
    :Q_min_inner_OuterControl,
    :Q_max_inner_OuterControl,
    :V_min_OuterControl,
    :V_max_OuterControl,
    :K_qp_OuterControl,
    :K_qi_OuterControl]

#REACTIVE POWER CONTROL
get_n_params(::PSY.ReactivePowerDroop) = 2
get_params(x::PSY.ReactivePowerDroop) = [PSY.get_kq(x), PSY.get_ωf(x)]
get_params_symbol(x::PSY.ReactivePowerDroop) = [:kq_OuterControl, :ωf_OuterControl]
#INNER CONTROL
get_n_params(::PSY.VoltageModeControl) = 10
get_params(x::PSY.VoltageModeControl) =
    [PSY.get_kpv(x), PSY.get_kiv(x), PSY.get_kffv(x), PSY.get_rv(x), PSY.get_lv(x),
        PSY.get_kpc(x), PSY.get_kic(x), PSY.get_kffi(x), PSY.get_ωad(x), PSY.get_kad(x)]
get_params_symbol(x::PSY.VoltageModeControl) =
    [
        :kpv_OuterControl,
        :kiv_OuterControl,
        :kffv_OuterControl,
        :rv_OuterControl,
        :lv_OuterControl,
        :kpc_OuterControl,
        :kic_OuterControl,
        :kffi_OuterControl,
        :ωad_OuterControl,
        :kad_OuterControl,
    ]

get_n_params(::PSY.RECurrentControlB) = 13
get_params(x::PSY.RECurrentControlB) =
    [
        PSY.get_Vdip_lim(x)[1],
        PSY.get_Vdip_lim(x)[2],
        PSY.get_T_rv(x),
        PSY.get_dbd_pnts(x)[1],
        PSY.get_dbd_pnts(x)[2],
        PSY.get_K_qv(x),
        PSY.get_Iqinj_lim(x)[1],
        PSY.get_Iqinj_lim(x)[2],
        PSY.get_V_ref0(x),
        PSY.get_K_vp(x),
        PSY.get_K_vi(x),
        PSY.get_T_iq(x),
        PSY.get_I_max(x)]
get_params_symbol(x::PSY.RECurrentControlB) =
    [:Vdip_min_OuterControl,
        :Vdip_max_OuterControl,
        :T_rv_OuterControl,
        :dbd_pnts_1_OuterControl,
        :dbd_pnts_2_OuterControl,
        :K_qv_OuterControl,
        :Iqinj_min_OuterControl,
        :Iqinj_max_OuterControl,
        :V_ref0_OuterControl,
        :K_vp_OuterControl,
        :K_vi_OuterControl,
        :T_iq_OuterControl,
        :I_max_OuterControl]

#DC SOURCE 
get_n_params(::PSY.FixedDCSource) = 1
get_params(x::PSY.FixedDCSource) = [PSY.get_voltage(x)]
get_params_symbol(x::PSY.FixedDCSource) = [:voltage_DCSource]
#FREQ ESTIMATOR
get_n_params(::PSY.KauraPLL) = 3
get_params(x::PSY.KauraPLL) = [PSY.get_ω_lp(x), PSY.get_kp_pll(x), PSY.get_ki_pll(x)]
get_params_symbol(::PSY.KauraPLL) =
    [:ω_lp_FrequencyEstimator, :kp_pll_FrequencyEstimator, :ki_pll_FrequencyEstimator]
get_n_params(::PSY.FixedFrequency) = 1
get_params(x::PSY.FixedFrequency) = [PSY.get_frequency(x)]
get_params_symbol(::PSY.FixedFrequency) = [:frequency_FrequencyEstimator]
#CONVERTER 
get_n_params(::PSY.AverageConverter) = 0
get_params(x::PSY.AverageConverter) = Float64[]
get_params_symbol(::PSY.AverageConverter) = Symbol[]
get_n_params(x::PSY.RenewableEnergyConverterTypeA) = 17
get_params(x::PSY.RenewableEnergyConverterTypeA) = [PSY.get_T_g(x),
    PSY.get_Rrpwr(x),
    PSY.get_Brkpt(x),
    PSY.get_Zerox(x),
    PSY.get_Lvpl1(x),
    PSY.get_Vo_lim(x),
    PSY.get_Lv_pnts(x)[1],
    PSY.get_Lv_pnts(x)[2],
    PSY.get_Io_lim(x),
    PSY.get_T_fltr(x),
    PSY.get_K_hv(x),
    PSY.get_Iqr_lims(x)[1],
    PSY.get_Iqr_lims(x)[2],
    PSY.get_Accel(x),
    PSY.get_Q_ref(x),
    PSY.get_R_source(x),
    PSY.get_X_source(x)]
get_params_symbol(::PSY.RenewableEnergyConverterTypeA) = [:T_g_Converter,
    :Rrpwr_Converter,
    :Brkpt_Converter,
    :Zerox_Converter,
    :Lvpl1_Converter,
    :Vo_lim_Converter,
    :Lv_pnt0_Converter,
    :Lv_pnt1_Converter,
    :Io_lim_Converter,
    :T_fltr_cnv_Converter,  #modified to make unique
    :K_hv_Converter,
    :Iqr_min_Converter,
    :Iqr_max_Converter,
    :Accel_Converter,
    :Q_ref_cnv_Converter,     #modified to make unique
    :R_source_Converter,
    :X_source_Converter]

#GENERATORS
get_n_params(g::PSY.DynamicGenerator) =
    3 +
    get_n_params(PSY.get_machine(g)) + get_n_params(PSY.get_shaft(g)) +
    get_n_params(PSY.get_avr(g)) + get_n_params(PSY.get_prime_mover(g)) +
    get_n_params(PSY.get_pss(g))
function get_params(g::PSY.DynamicGenerator)
    refs = [
        PSY.get_V_ref(g),
        PSY.get_ω_ref(g),
        PSY.get_P_ref(g),
    ]
    vcat(
        refs,
        get_params(PSY.get_machine(g)),
        get_params(PSY.get_shaft(g)),
        get_params(PSY.get_avr(g)),
        get_params(PSY.get_prime_mover(g)),
        get_params(PSY.get_pss(g)),
    )
end
get_params_symbol(g::PSY.DynamicGenerator) = vcat(
    [:V_ref, :ω_ref, :P_ref],
    get_params_symbol(PSY.get_machine(g)),
    get_params_symbol(PSY.get_shaft(g)),
    get_params_symbol(PSY.get_avr(g)),
    get_params_symbol(PSY.get_prime_mover(g)),
    get_params_symbol(PSY.get_pss(g)),
)

#MACHINES 
get_n_params(::PSY.BaseMachine) = 3
get_params(x::PSY.BaseMachine) = [PSY.get_R(x), PSY.get_Xd_p(x), PSY.get_eq_p(x)]
get_params_symbol(::PSY.BaseMachine) = [:R_Machine, :Xd_p_Machine, :eq_p_Machine]
get_n_params(::PSY.OneDOneQMachine) = 7
get_params(x::PSY.OneDOneQMachine) = [
    PSY.get_R(x),
    PSY.get_Xd(x),
    PSY.get_Xq(x),
    PSY.get_Xd_p(x),
    PSY.get_Xq_p(x),
    PSY.get_Td0_p(x),
    PSY.get_Tq0_p(x),
]
get_params_symbol(x::PSY.OneDOneQMachine) = [
    :R_Machine,
    :Xd_Machine,
    :Xq_Machine,
    :Xd_p_Machine,
    :Xq_p_Machine,
    :Td0_p_Machine,
    :Tq0_p_Machine,
]
get_n_params(::PSY.MarconatoMachine) = 14
get_params(x::PSY.MarconatoMachine) =
    [PSY.get_R(x), PSY.get_Xd(x), PSY.get_Xq(x), PSY.get_Xd_p(x), PSY.get_Xq_p(x),
        PSY.get_Xd_pp(x), PSY.get_Xq_pp(x),
        PSY.get_Td0_p(x), PSY.get_Tq0_p(x), PSY.get_Td0_pp(x), PSY.get_Tq0_pp(x),
        PSY.get_T_AA(x), PSY.get_γd(x), PSY.get_γq(x)]
get_params_symbol(x::PSY.MarconatoMachine) = [
    :R_Machine,
    :Xd_Machine,
    :Xq_Machine,
    :Xd_p_Machine,
    :Xq_p_Machine,
    :Xd_pp_Machine,
    :Xq_pp_Machine,
    :Td0_p_Machine,
    :Tq0_p_Machine,
    :Td0_pp_Machine,
    :Tq0_pp_Machine,
    :T_AA_Machine,
    :γd_Machine,
    :γq_Machine,
]
get_n_params(::PSY.AndersonFouadMachine) = 11
get_params(x::PSY.AndersonFouadMachine) = [
    PSY.get_R(x),
    PSY.get_Xd(x),
    PSY.get_Xq(x),
    PSY.get_Xd_p(x),
    PSY.get_Xq_p(x),
    PSY.get_Xd_pp(x),
    PSY.get_Xq_pp(x),
    PSY.get_Td0_p(x),
    PSY.get_Tq0_p(x),
    PSY.get_Td0_pp(x),
    PSY.get_Tq0_pp(x)]

get_params_symbol(x::PSY.AndersonFouadMachine) = [:R_Machine,
    :Xd_Machine,
    :Xq_Machine,
    :Xd_p_Machine,
    :Xq_p_Machine,
    :Xd_pp_Machine,
    :Xq_pp_Machine,
    :Td0_p_Machine,
    :Tq0_p_Machine,
    :Td0_pp_Machine,
    :Tq0_pp_Machine]

#NOTE: Saturation not considered as paramters
get_n_params(
    x::Union{PSY.RoundRotorMachine, PSY.RoundRotorExponential, PSY.RoundRotorQuadratic},
) = 16
get_params(
    x::Union{PSY.RoundRotorMachine, PSY.RoundRotorExponential, PSY.RoundRotorQuadratic},
) = [
    PSY.get_R(x),
    PSY.get_Td0_p(x),
    PSY.get_Td0_pp(x),
    PSY.get_Tq0_p(x),
    PSY.get_Tq0_pp(x),
    PSY.get_Xd(x),
    PSY.get_Xq(x),
    PSY.get_Xd_p(x),
    PSY.get_Xq_p(x),
    PSY.get_Xd_pp(x),
    PSY.get_Xl(x),
    PSY.get_γ_d1(x),
    PSY.get_γ_q1(x),
    PSY.get_γ_d2(x),
    PSY.get_γ_q2(x),
    PSY.get_γ_qd(x),
]
get_params_symbol(
    ::Union{PSY.RoundRotorMachine, PSY.RoundRotorExponential, PSY.RoundRotorQuadratic},
) = [
    :R_Machine,
    :Td0_p_Machine,
    :Td0_pp_Machine,
    :Tq0_p_Machine,
    :Tq0_pp_Machine,
    :Xd_Machine,
    :Xq_Machine,
    :Xd_p_Machine,
    :Xq_p_Machine,
    :Xd_pp_Machine,
    :Xl_Machine,
    :γ_d1_Machine,
    :γ_q1_Machine,
    :γ_d2_Machine,
    :γ_q2_Machine,
    :γ_qd_Machine,
]
get_n_params(
    ::Union{PSY.SalientPoleMachine, PSY.SalientPoleExponential, PSY.SalientPoleQuadratic},
) = 12
get_params(
    x::Union{PSY.SalientPoleMachine, PSY.SalientPoleExponential, PSY.SalientPoleQuadratic},
) = [
    PSY.get_R(x),
    PSY.get_Td0_p(x),
    PSY.get_Td0_pp(x),
    PSY.get_Tq0_pp(x),
    PSY.get_Xd(x),
    PSY.get_Xq(x),
    PSY.get_Xd_p(x),
    PSY.get_Xd_pp(x),
    PSY.get_Xl(x),
    PSY.get_γ_d1(x),
    PSY.get_γ_q1(x),
    PSY.get_γ_d2(x),
]
get_params_symbol(
    ::Union{PSY.SalientPoleMachine, PSY.SalientPoleExponential, PSY.SalientPoleQuadratic},
) = [
    :R_Machine,
    :Td0_p_Machine,
    :Td0_pp_Machine,
    :Tq0_pp_Machine,
    :Xd_Machine,
    :Xq_Machine,
    :Xd_p_Machine,
    :Xd_pp_Machine,
    :Xl_Machine,
    :γ_d1_Machine,
    :γ_q1_Machine,
    :γ_d2_Machine,
]

#SHAFTS
get_n_params(::PSY.SingleMass) = 2
get_params(x::PSY.SingleMass) = [PSY.get_H(x), PSY.get_D(x)]
get_params_symbol(::PSY.SingleMass) = [:H_Shaft, :D_Shaft]
get_n_params(::PSY.FiveMassShaft) = 18
get_params(x::PSY.FiveMassShaft) = [
    PSY.get_H(x),
    PSY.get_H_hp(x),
    PSY.get_H_ip(x),
    PSY.get_H_lp(x),
    PSY.get_H_ex(x),
    PSY.get_D(x),
    PSY.get_D_hp(x),
    PSY.get_D_ip(x),
    PSY.get_D_lp(x),
    PSY.get_D_ex(x),
    PSY.get_D_12(x),
    PSY.get_D_23(x),
    PSY.get_D_34(x),
    PSY.get_D_45(x),
    PSY.get_K_hp(x),
    PSY.get_K_ip(x),
    PSY.get_K_lp(x),
    PSY.get_K_ex(x),
]

get_params_symbol(::PSY.FiveMassShaft) = [
    :H_Shaft,
    :H_hp_Shaft,
    :H_ip_Shaft,
    :H_lp_Shaft,
    :H_ex_Shaft,
    :D_Shaft,
    :D_hp_Shaft,
    :D_ip_Shaft,
    :D_lp_Shaft,
    :D_ex_Shaft,
    :D_12_Shaft,
    :D_23_Shaft,
    :D_34_Shaft,
    :D_45_Shaft,
    :K_hp_Shaft,
    :K_ip_Shaft,
    :K_lp_Shaft,
    :K_ex_Shaft]

#AVRS 
get_n_params(::PSY.AVRFixed) = 0
get_params(x::PSY.AVRFixed) = Float64[]
get_params_symbol(::PSY.AVRFixed) = Symbol[]
get_n_params(::PSY.AVRSimple) = 1
get_params(x::PSY.AVRSimple) = [PSY.get_Kv(x)]
get_params_symbol(::PSY.AVRSimple) = [:Kv_AVR]
get_n_params(::PSY.AVRTypeI) = 9
get_params(x::PSY.AVRTypeI) = [
    PSY.get_Ka(x),
    PSY.get_Ke(x),
    PSY.get_Kf(x),
    PSY.get_Ta(x),
    PSY.get_Te(x),
    PSY.get_Tf(x),
    PSY.get_Tr(x),
    PSY.get_Ae(x),
    PSY.get_Be(x),
]
get_params_symbol(::PSY.AVRTypeI) =
    [:Ka_AVR, :Ke_AVR, :Kf_AVR, :Ta_AVR, :Te_AVR, :Tf_AVR, :Tr_AVR, :Ae_AVR, :Be_AVR]
get_n_params(::PSY.SEXS) = 6
get_params(x::PSY.SEXS) = [
    PSY.get_Ta_Tb(x),
    PSY.get_Tb(x),
    PSY.get_K(x),
    PSY.get_Te(x),
    PSY.get_V_lim(x)[1],
    PSY.get_V_lim(x)[2],
]
get_params_symbol(::PSY.SEXS) =
    Symbol[:Ta_Tb_AVR, :Tb_AVR, :K_AVR, :Te_AVR, :V_min_AVR, :V_max_AVR]
get_n_params(::PSY.AVRTypeII) = 11
get_params(x::PSY.AVRTypeII) =
    [PSY.get_K0(x),
        PSY.get_T1(x),
        PSY.get_T2(x),
        PSY.get_T3(x),
        PSY.get_T4(x),
        PSY.get_Te(x),
        PSY.get_Tr(x),
        PSY.get_Va_lim(x)[1],
        PSY.get_Va_lim(x)[2],
        PSY.get_Ae(x),
        PSY.get_Be(x)]
get_params_symbol(::PSY.AVRTypeII) =
    [:K0_AVR,
        :T1_AVR,
        :T2_AVR,
        :T3_AVR,
        :T4_AVR,
        :Te_AVR,
        :Tr_AVR,
        :Va_min_AVR,
        :Va_max_AVR,
        :Ae_AVR,
        :Be_AVR]
get_n_params(::PSY.ESAC1A) = 15
get_params(x::PSY.ESAC1A) = [
    PSY.get_Tr(x),
    PSY.get_Tb(x),
    PSY.get_Tc(x),
    PSY.get_Ka(x),
    PSY.get_Ta(x),
    PSY.get_Va_lim(x)[1],
    PSY.get_Va_lim(x)[2],
    PSY.get_Te(x),
    PSY.get_Kf(x),
    PSY.get_Tf(x),
    PSY.get_Kc(x),
    PSY.get_Kd(x),
    PSY.get_Ke(x),
    PSY.get_Vr_lim(x)[1],
    PSY.get_Vr_lim(x)[2]]
get_params_symbol(::PSY.ESAC1A) = [:Tr,
    :Tb_AVR,
    :Tc_AVR,
    :Ka_AVR,
    :Ta_AVR,
    :Va_min_AVR,
    :Va_max_AVR,
    :Te_AVR,
    :Kf_AVR,
    :Tf_AVR,
    :Kc_AVR,
    :Kd_AVR,
    :Ke_AVR,
    :Vr_min_AVR,
    :Vr_max_AVR]
get_n_params(::PSY.EXST1) = 12
get_params(x::PSY.EXST1) = [PSY.get_Tr(x),
    PSY.get_Vi_lim(x)[1],
    PSY.get_Vi_lim(x)[2],
    PSY.get_Tc(x),
    PSY.get_Tb(x),
    PSY.get_Ka(x),
    PSY.get_Ta(x),
    PSY.get_Vr_lim(x)[1],
    PSY.get_Vr_lim(x)[2],
    PSY.get_Kc(x),
    PSY.get_Kf(x),
    PSY.get_Tf(x)]
get_params_symbol(::PSY.EXST1) = [
    :Tr_AVR,
    :Vi_min_AVR,
    :Vi_max_AVR,
    :Tc_AVR,
    :Tb_AVR,
    :Ka_AVR,
    :Ta_AVR,
    :Vr_min_AVR,
    :Vr_max_AVR,
    :Kc_AVR,
    :Kf_AVR,
    :Tf_AVR,
]

get_n_params(::PSY.EXAC1) = 13
get_params(x::PSY.EXAC1) = [
    PSY.get_Tr(x),
    PSY.get_Tb(x),
    PSY.get_Tc(x),
    PSY.get_Ka(x),
    PSY.get_Ta(x),
    PSY.get_Vr_lim(x)[1],
    PSY.get_Vr_lim(x)[2],
    PSY.get_Te(x),
    PSY.get_Kf(x),
    PSY.get_Tf(x),
    PSY.get_Kc(x),
    PSY.get_Kd(x),
    PSY.get_Ke(x)]

get_params_symbol(::PSY.EXAC1) = [
    :Tr_AVR,
    :Tb_AVR,
    :Tc_AVR,
    :Ka_AVR,
    :Ta_AVR,
    :Vr_min_AVR,
    :Vr_max_AVR,
    :Te_AVR,
    :Kf_AVR,
    :Tf_AVR,
    :Kc_AVR,
    :Kd_AVR,
    :Ke_AVR]

#TurbineGov
get_n_params(::PSY.TGFixed) = 1
get_params(x::PSY.TGFixed) = [PSY.get_efficiency(x)]
get_params_symbol(::PSY.TGFixed) = [:efficiency_TurbineGov]
get_n_params(::PSY.TGTypeII) = 3
get_params(x::PSY.TGTypeII) = [PSY.get_R(x), PSY.get_T1(x), PSY.get_T2(x)]
get_params_symbol(::PSY.TGTypeII) = [:R_tg_TurbineGovR, :T1_TurbineGov, :T2_TurbineGov]
get_n_params(::PSY.GasTG) = 9
get_params(x::PSY.GasTG) = [
    PSY.get_R(x),
    PSY.get_T1(x),
    PSY.get_T2(x),
    PSY.get_T3(x),
    PSY.get_AT(x),
    PSY.get_Kt(x),
    PSY.get_V_lim(x)[1],
    PSY.get_V_lim(x)[2],
    PSY.get_D_turb(x),
]
get_params_symbol(::PSY.GasTG) = [
    :R_tg_TurbineGov,
    :T1_TurbineGov,
    :T2_TurbineGov,
    :T3_TurbineGov,
    :AT_TurbineGov,
    :Kt_TurbineGov,
    :V_min_TurbineGov,
    :V_max_TurbineGov,
    :D_turb_TurbineGov,
]
get_n_params(::PSY.TGTypeI) = 8
get_params(x::PSY.TGTypeI) = [PSY.get_R(x),
    PSY.get_Ts(x),
    PSY.get_Tc(x),
    PSY.get_T3(x),
    PSY.get_T4(x),
    PSY.get_T5(x),
    PSY.get_valve_position_limits(x)[1],
    PSY.get_valve_position_limits(x)[2],
]

get_params_symbol(::PSY.TGTypeI) = [:R_tg_TurbineGov,
    :Ts_TurbineGov,
    :Tc_TurbineGov,
    :T3_TurbineGov,
    :T4_TurbineGov,
    :T5_TurbineGov,
    :valve_position_min_TurbineGov,
    :valve_position_max_TurbineGov,
]
get_n_params(::PSY.SteamTurbineGov1) = 7
get_params(x::PSY.SteamTurbineGov1) = [PSY.get_R(x),
    PSY.get_T1(x),
    PSY.get_valve_position_limits(x)[1],
    PSY.get_valve_position_limits(x)[2],
    PSY.get_T2(x),
    PSY.get_T3(x),
    PSY.get_D_T(x)]
get_params_symbol(::PSY.SteamTurbineGov1) = [:R_tg,
    :T1_TurbineGov,
    :valve_position_min_TurbineGov,
    :valve_position_max_TurbineGov,
    :T2_TurbineGov,
    :T3_TurbineGov,
    :D_T_TurbineGov]
get_n_params(::PSY.HydroTurbineGov) = 12
get_params(x::PSY.HydroTurbineGov) = [PSY.get_R(x),
    PSY.get_r(x),
    PSY.get_Tr(x),
    PSY.get_Tf(x),
    PSY.get_Tg(x),
    PSY.get_VELM(x),
    PSY.get_gate_position_limits(x)[1],
    PSY.get_gate_position_limits(x)[2],
    PSY.get_Tw(x),
    PSY.get_At(x),
    PSY.get_D_T(x),
    PSY.get_q_nl(x)]
get_params_symbol(::PSY.HydroTurbineGov) = [
    :R_TurbineGov,
    :r_TurbineGov,
    :Tr_TurbineGov,
    :Tf_TurbineGov,
    :Tg_TurbineGov,
    :VELM_TurbineGov,
    :G_min_TurbineGov,
    :G_max_TurbineGov,
    :Tw_TurbineGov,
    :At_TurbineGov,
    :D_T_TurbineGov,
    :q_nl_TurbineGov]

#PSS
get_n_params(::PSY.PSSFixed) = 1
get_params(x::PSY.PSSFixed) = [PSY.get_V_pss(x)]
get_params_symbol(::PSY.PSSFixed) = [:V_pss_PSS]

get_n_params(::PSY.STAB1) = 7
get_params(x::PSY.STAB1) = [PSY.get_KT(x),
    PSY.get_T(x),
    PSY.get_T1T3(x),
    PSY.get_T3(x),
    PSY.get_T2T4(x),
    PSY.get_T4(x),
    PSY.get_H_lim(x)]
get_params_symbol(::PSY.STAB1) =
    [:KT_PSS,
        :T_PSS,
        :T1T3_PSS,
        :T3_PSS,
        :T2T4_PSS,
        :T4_PSS,
        :H_lim_PSS]

#STATIC INJECTION
get_n_params(::PSY.StaticInjection) = 1
get_params(x::PSY.StaticInjection) = [PSY.get_reactive_power(x)]
get_params_symbol(::PSY.StaticInjection) = [:Q_ref]

get_n_params(::PSY.StandardLoad) = 1
get_params(x::PSY.StandardLoad) = [PF.get_total_q(x)]
get_params_symbol(::PSY.StandardLoad) = [:Q_ref]
#SOURCE 
get_n_params(::PSY.Source) = 4
get_params(x::PSY.Source) = [
    PSY.get_R_th(x),
    PSY.get_X_th(x),
    PSY.get_internal_voltage(x),
    PSY.get_internal_angle(x),
]
#Parameters not implemented for PVS - requires change in PSY Struct to have information required to construct and deconstruct parameter vector
get_n_params(x::PSY.PeriodicVariableSource) = 4
get_params(x::PSY.PeriodicVariableSource) = Float64[0.0, 0.0, 0.0, 0.0]

#DYNAMIC LOADS 
get_n_params(::PSY.ActiveConstantPowerLoad) = 18
get_params(x::PSY.ActiveConstantPowerLoad) =
    [
        0.0,
        0.0,
        0.0,
        0.0,
        PSY.get_r_load(x),
        PSY.get_c_dc(x),
        PSY.get_rf(x),
        PSY.get_lf(x),
        PSY.get_cf(x),
        PSY.get_rg(x),
        PSY.get_lg(x),
        PSY.get_kp_pll(x),
        PSY.get_ki_pll(x),
        PSY.get_kpv(x),
        PSY.get_kiv(x),
        PSY.get_kpc(x),
        PSY.get_kic(x),
        PSY.get_base_power(x)]
get_params_symbol(::PSY.ActiveConstantPowerLoad) =
    [
        :Q_ref,
        :V_ref,
        :ω_ref,
        :P_ref,
        :r_load,
        :c_dc,
        :rf,
        :lf,
        :cf,
        :rg,
        :lg,
        :kp_pll,
        :ki_pll,
        :kpv,
        :kiv,
        :kpc,
        :kic,
        :base_power]

get_n_params(::PSY.SingleCageInductionMachine) = 18
get_params(x::PSY.SingleCageInductionMachine) = [
    0.0,#refs unused
    0.0,#refs unused
    0.0, #refs unused
    0.0,  #refs unused
    PSY.get_R_s(x),
    PSY.get_R_r(x),
    PSY.get_X_ls(x),
    PSY.get_X_lr(x),
    PSY.get_X_m(x),
    PSY.get_H(x),
    PSY.get_A(x),
    PSY.get_B(x),
    PSY.get_base_power(x),
    PSY.get_C(x),
    PSY.get_τ_ref(x),
    PSY.get_B_shunt(x),
    PSY.get_X_ad(x),
    PSY.get_X_aq(x)]
get_params_symbol(::PSY.SingleCageInductionMachine) = [
    :Q_ref,
    :V_ref,
    :ω_ref,
    :P_ref,
    :R_s,
    :R_r,
    :X_ls,
    :X_lr,
    :X_m,
    :H,
    :A,
    :B,
    :base_power,
    :C,
    :τ_ref,
    :B_shunt,
    :X_ad,
    :X_aq]
get_n_params(::PSY.SimplifiedSingleCageInductionMachine) = 19
get_params(x::PSY.SimplifiedSingleCageInductionMachine) = [
    0.0,#refs unused
    0.0,#refs unused
    0.0, #refs unused
    0.0,  #refs unused
    PSY.get_R_s(x),
    PSY.get_R_r(x),
    PSY.get_X_ls(x),
    PSY.get_X_lr(x),
    PSY.get_X_m(x),
    PSY.get_H(x),
    PSY.get_A(x),
    PSY.get_B(x),
    PSY.get_base_power(x),
    PSY.get_C(x),
    PSY.get_τ_ref(x),
    PSY.get_B_shunt(x),
    PSY.get_X_ss(x),
    PSY.get_X_rr(x),
    PSY.get_X_p(x)]
get_params_symbol(::PSY.SimplifiedSingleCageInductionMachine) = [
    :Q_ref,
    :V_ref,
    :ω_ref,
    :P_ref,
    :R_s,
    :R_r,
    :X_ls,
    :X_lr,
    :X_m,
    :H,
    :A,
    :B,
    :base_power,
    :C,
    :τ_ref,
    :B_shunt,
    :X_ss,
    :X_rr,
    :X_p]

get_n_params(x::PSY.AggregateDistributedGenerationA) = 33
get_params(x::PSY.AggregateDistributedGenerationA) = [PSY.get_Q_ref(x),
    PSY.get_V_ref(x),
    PSY.get_ω_ref(x),
    PSY.get_P_ref(x),
    PSY.get_T_rv(x),
    PSY.get_Trf(x),
    PSY.get_dbd_pnts(x)[1],
    PSY.get_dbd_pnts(x)[2],
    PSY.get_K_qv(x),
    PSY.get_Tp(x),
    PSY.get_T_iq(x),
    PSY.get_D_dn(x),
    PSY.get_D_up(x),
    PSY.get_fdbd_pnts(x)[1],
    PSY.get_fdbd_pnts(x)[2],
    PSY.get_fe_lim(x)[1],
    PSY.get_fe_lim(x)[2],
    PSY.get_P_lim(x)[1],
    PSY.get_P_lim(x)[2],
    PSY.get_dP_lim(x)[1],
    PSY.get_dP_lim(x)[2],
    PSY.get_Tpord(x),
    PSY.get_Kpg(x),
    PSY.get_Kig(x),
    PSY.get_I_max(x),
    PSY.get_Tg(x),
    PSY.get_rrpwr(x),
    PSY.get_Tv(x),
    PSY.get_Vpr(x),
    PSY.get_Iq_lim(x)[1],
    PSY.get_Iq_lim(x)[2],
    PSY.get_base_power(x),
    PSY.get_Pfa_ref(x)]

get_n_params(x::PSY.CSVGN1) = 17
get_params(x::PSY.CSVGN1) = [
    0.0,
    0.0,
    0.0,
    0.0,
    PSY.get_K(x),
    PSY.get_T1(x),
    PSY.get_T2(x),
    PSY.get_T3(x),
    PSY.get_T4(x),
    PSY.get_T5(x),
    PSY.get_Rmin(x),
    PSY.get_Vmax(x),
    PSY.get_Vmin(x),
    PSY.get_CBase(x),
    PSY.get_base_power(x),
    PSY.get_R_th(x),
    PSY.get_X_th(x)]
