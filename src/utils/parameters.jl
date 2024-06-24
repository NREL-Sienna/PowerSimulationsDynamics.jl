
@enum PARAM_TYPES begin
    DEVICE_PARAM = 0
    DEVICE_SETPOINT = 1
    NETWORK_PARAM = 2
    NETWORK_SETPOINT = 3
end

struct ParamsMetadata
    type::PARAM_TYPES
    in_mass_matrix::Bool
    impacts_ic::Bool
end
Base.length(::ParamsMetadata) = 1

function get_params(x)
    @error "Parameters not yet defined for component: $(typeof(x))"
    (;)
end

get_params_metadata(::T) where {T <: PSY.DynamicComponent} = (;)
get_params(::PSY.ActivePowerControl) = (;)
get_params_metadata(::PSY.ActivePowerControl) = (;)
get_params(::PSY.ReactivePowerControl) = (;)
get_params_metadata(::PSY.ReactivePowerControl) = (;)

get_params(
    d::DynamicWrapper{T},
) where {T <: Union{PSY.DynamicGenerator, PSY.DynamicInverter}} = get_params(get_device(d))
get_params_metadata(
    d::DynamicWrapper{T},
) where {T <: Union{PSY.DynamicGenerator, PSY.DynamicInverter}} =
    get_params_metadata(get_device(d))

get_params(d::DynamicWrapper) = get_params(get_device(d))
get_params_metadata(d::DynamicWrapper) = get_params_metadata(get_device(d))

get_params(d::StaticWrapper) = get_params(get_device(d))
get_params_metadata(d::StaticWrapper) = get_params_metadata(get_device(d))

get_params(x::BranchWrapper) = get_params(get_branch(x))
get_params_metadata(x::BranchWrapper) = get_params_metadata(get_branch(x))

get_params(x::PSY.DynamicBranch) = get_params(PSY.get_branch(x))
get_params_metadata(x::PSY.DynamicBranch) = get_params_metadata(PSY.get_branch(x))

get_params(x::PSY.Line) = (r = PSY.get_r(x), x = PSY.get_x(x))
get_params_metadata(::PSY.Line) = (
    r = ParamsMetadata(NETWORK_PARAM, false, true),
    x = ParamsMetadata(NETWORK_PARAM, false, true),
)
get_params(::StaticLoadWrapper) = (;)
get_params_metadata(::StaticLoadWrapper) = (;)
########### INVERTERS #############
function get_params(g::PSY.DynamicInverter)
    (
        Converter = get_params(PSY.get_converter(g)),
        DCSource = get_params(PSY.get_dc_source(g)),
        Filter = get_params(PSY.get_filter(g)),
        FrequencyEstimator = get_params(PSY.get_freq_estimator(g)),
        InnerControl = get_params(PSY.get_inner_control(g)),
        OuterControl = get_params(PSY.get_outer_control(g)),
    )
end
function get_params_metadata(g::PSY.DynamicInverter)
    (
        Converter = get_params_metadata(PSY.get_converter(g)),
        DCSource = get_params_metadata(PSY.get_dc_source(g)),
        Filter = get_params_metadata(PSY.get_filter(g)),
        FrequencyEstimator = get_params_metadata(PSY.get_freq_estimator(g)),
        InnerControl = get_params_metadata(PSY.get_inner_control(g)),
        OuterControl = get_params_metadata(PSY.get_outer_control(g)),
    )
end

#FILTERS 
get_params(x::PSY.LCLFilter) =
    (
        lf = PSY.get_lf(x),
        rf = PSY.get_rf(x),
        cf = PSY.get_cf(x),
        lg = PSY.get_lg(x),
        rg = PSY.get_rg(x),
    )
get_params_metadata(::PSY.LCLFilter) = (
    lf = ParamsMetadata(DEVICE_PARAM, true, true),
    rf = ParamsMetadata(DEVICE_PARAM, false, true),
    cf = ParamsMetadata(DEVICE_PARAM, true, true),
    lg = ParamsMetadata(DEVICE_PARAM, true, true),
    rg = ParamsMetadata(DEVICE_PARAM, false, true),
)
get_params(x::PSY.RLFilter) = (rf = PSY.get_rf(x), lf = PSY.get_lf(x))
get_params_metadata(::PSY.RLFilter) = (
    rf = ParamsMetadata(DEVICE_PARAM, false, true),
    lf = ParamsMetadata(DEVICE_PARAM, false, true),
)

#OUTER CONTROL,                                                                                                                                                                
get_params(x::PSY.OuterControl) = (
    ActivePowerControl = get_params(PSY.get_active_power_control(x)),
    ReactivePowerControl = get_params(PSY.get_reactive_power_control(x)),
)
get_params_metadata(x::PSY.OuterControl) = (
    ActivePowerControl = get_params_metadata(PSY.get_active_power_control(x)),
    ReactivePowerControl = get_params_metadata(PSY.get_reactive_power_control(x)),
)
#ACTIVE POWER CONTROL
get_params(x::PSY.VirtualInertia) =
    (Ta = PSY.get_Ta(x), kd = PSY.get_kd(x), kω = PSY.get_kω(x))
get_params_metadata(::PSY.VirtualInertia) = (
    Ta = ParamsMetadata(DEVICE_PARAM, false, false),
    kd = ParamsMetadata(DEVICE_PARAM, false, false),
    kω = ParamsMetadata(DEVICE_PARAM, false, false),
)
#Note: Removed fbdd_pnts from parameters because it is not a NamedTuple
get_params(x::PSY.ActiveRenewableControllerAB) = (
    K_pg = PSY.get_K_pg(x),
    K_ig = PSY.get_K_ig(x),
    T_p = PSY.get_T_p(x),
    fe_lim = PSY.get_fe_lim(x),
    P_lim = PSY.get_P_lim(x),
    T_g = PSY.get_T_g(x),
    D_dn = PSY.get_D_dn(x),
    D_up = PSY.get_D_up(x),
    dP_lim = PSY.get_dP_lim(x),
    P_lim_inner = PSY.get_P_lim_inner(x),
    T_pord = PSY.get_T_pord(x),
)
get_params_metadata(::PSY.ActiveRenewableControllerAB) = (
    K_pg = ParamsMetadata(DEVICE_PARAM, false, false),
    K_ig = ParamsMetadata(DEVICE_PARAM, false, true),
    T_p = ParamsMetadata(DEVICE_PARAM, false, true),
    fe_lim = (
        min = ParamsMetadata(DEVICE_PARAM, false, false),
        max = ParamsMetadata(DEVICE_PARAM, false, false),
    ),
    P_lim = (
        min = ParamsMetadata(DEVICE_PARAM, false, false),
        max = ParamsMetadata(DEVICE_PARAM, false, false),
    ),
    T_g = ParamsMetadata(DEVICE_PARAM, false, false),
    D_dn = ParamsMetadata(DEVICE_PARAM, false, false),
    D_up = ParamsMetadata(DEVICE_PARAM, false, false),
    dP_lim = (
        min = ParamsMetadata(DEVICE_PARAM, false, false),
        max = ParamsMetadata(DEVICE_PARAM, false, false),
    ),
    P_lim_innem_inner = (
        min = ParamsMetadata(DEVICE_PARAM, false, false),
        max = ParamsMetadata(DEVICE_PARAM, false, false),
    ),
    T_pord = ParamsMetadata(DEVICE_PARAM, false, false),
)
#REACTIVE POWER CONTROL
get_params(x::PSY.ReactivePowerDroop) = (kq = PSY.get_kq(x), ωf = PSY.get_ωf(x))
get_params_metadata(::PSY.ReactivePowerDroop) = (
    kq = ParamsMetadata(DEVICE_PARAM, false, false),
    ωf = ParamsMetadata(DEVICE_PARAM, false, false),
)
get_params(x::PSY.ReactiveRenewableControllerAB) = (
    T_fltr = PSY.get_T_fltr(x),
    K_p = PSY.get_K_p(x),
    K_i = PSY.get_K_i(x),
    T_ft = PSY.get_T_ft(x),
    T_fv = PSY.get_T_fv(x),
    V_frz = PSY.get_V_frz(x),
    R_c = PSY.get_R_c(x),
    X_c = PSY.get_X_c(x),
    K_c = PSY.get_K_c(x),
    e_lim = (min = PSY.get_e_lim(x)[:min], max = PSY.get_e_lim(x)[:max]),
    Q_lim = (min = PSY.get_Q_lim(x)[:min], max = PSY.get_Q_lim(x)[:max]),
    T_p = PSY.get_T_p(x),
    Q_lim_inner = (min = PSY.get_Q_lim_inner(x)[:min], max = PSY.get_Q_lim_inner(x)[:max]),
    V_lim = (min = PSY.get_V_lim(x)[:min], max = PSY.get_V_lim(x)[:max]),
    K_qp = PSY.get_K_qp(x),
    K_qi = PSY.get_K_qi(x),
)
get_params_metadata(::PSY.ReactiveRenewableControllerAB) =
    (T_fltr = ParamsMetadata(DEVICE_PARAM, false, false),
        K_p = ParamsMetadata(DEVICE_PARAM, false, false),
        K_i = ParamsMetadata(DEVICE_PARAM, false, true),
        T_ft = ParamsMetadata(DEVICE_PARAM, false, false),
        T_fv = ParamsMetadata(DEVICE_PARAM, false, false),
        V_frz = ParamsMetadata(DEVICE_PARAM, false, false),
        R_c = ParamsMetadata(DEVICE_PARAM, false, true),
        X_c = ParamsMetadata(DEVICE_PARAM, false, true),
        K_c = ParamsMetadata(DEVICE_PARAM, false, true),
        e_lim = (
            min = ParamsMetadata(DEVICE_PARAM, false, false),
            max = ParamsMetadata(DEVICE_PARAM, false, false),
        ),
        Q_lim = (
            min = ParamsMetadata(DEVICE_PARAM, false, false),
            max = ParamsMetadata(DEVICE_PARAM, false, false),
        ),
        T_p = ParamsMetadata(DEVICE_PARAM, false, true),
        Q_lim_inner = (
            min = ParamsMetadata(DEVICE_PARAM, false, false),
            max = ParamsMetadata(DEVICE_PARAM, false, false),
        ),
        V_lim = (
            min = ParamsMetadata(DEVICE_PARAM, false, false),
            max = ParamsMetadata(DEVICE_PARAM, false, false),
        ),
        K_qp = ParamsMetadata(DEVICE_PARAM, false, false),
        K_qi = ParamsMetadata(DEVICE_PARAM, false, true),
    )
#= 
= PSY.get_dbd_pnts(x)[1],
= PSY.get_dbd_pnts(x)[2],

ParamsMetadata(:dbd_pnts1_ReactivePowerControl, false, false, false, false),
ParamsMetadata(:dbd_pnts2_ReactivePowerControl, false, false, false, false),
 =#

#INNER CONTROL
get_params(x::PSY.VoltageModeControl) = (
    kpv = PSY.get_kpv(x),
    kiv = PSY.get_kiv(x),
    kffv = PSY.get_kffv(x),
    rv = PSY.get_rv(x),
    lv = PSY.get_lv(x),
    kpc = PSY.get_kpc(x),
    kic = PSY.get_kic(x),
    kffi = PSY.get_kffi(x),
    ωad = PSY.get_ωad(x),
    kad = PSY.get_kad(x),
)
get_params_metadata(::PSY.VoltageModeControl) = (
    kpv = ParamsMetadata(DEVICE_PARAM, false, true),
    kiv = ParamsMetadata(DEVICE_PARAM, false, true),
    kffv = ParamsMetadata(DEVICE_PARAM, false, true),
    rv = ParamsMetadata(DEVICE_PARAM, false, true),
    lv = ParamsMetadata(DEVICE_PARAM, false, true),
    kpc = ParamsMetadata(DEVICE_PARAM, false, true),
    kic = ParamsMetadata(DEVICE_PARAM, false, true),
    kffi = ParamsMetadata(DEVICE_PARAM, false, true),
    ωad = ParamsMetadata(DEVICE_PARAM, false, false),
    kad = ParamsMetadata(DEVICE_PARAM, false, true),
)
#= 
get_params(x::PSY.RECurrentControlB) = [
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
    PSY.get_I_max(x),
]
get_params_metadata(::PSY.RECurrentControlB) = [
    ParamsMetadata(:Vdip_min_InnerControl, false, false, false, false),
    ParamsMetadata(:Vdip_max_InnerControl, false, false, false, false),
    ParamsMetadata(:T_rv_InnerControl, true, false, false, false),
    ParamsMetadata(:dbd_pnts_1_InnerControl, false, false, false, false),
    ParamsMetadata(:dbd_pnts_2_InnerControl, false, false, false, false),
    ParamsMetadata(:K_qv_InnerControl, false, false, false, false),
    ParamsMetadata(:Iqinj_min_InnerControl, false, false, false, false),
    ParamsMetadata(:Iqinj_max_InnerControl, false, false, false, false),
    ParamsMetadata(:V_ref0_InnerControl, false, false, true, false),
    ParamsMetadata(:K_vp_InnerControl, false, false, false, false),
    ParamsMetadata(:K_vi_InnerControl, false, false, true, false),
    ParamsMetadata(:T_iq_InnerControl, true, false, false, false),
    ParamsMetadata(:I_max_InnerControl, false, false, false, false),
]
 =#
#DC SOURCE 
get_params(x::PSY.FixedDCSource) = (voltage = PSY.get_voltage(x),)
get_params_metadata(::PSY.FixedDCSource) =
    (voltage = ParamMetadata(DEVICE_PARAM, false, false),)

#FREQ ESTIMATOR
get_params(x::PSY.KauraPLL) = (
    ω_lp = PSY.get_ω_lp(x),
    kp_pll = PSY.get_kp_pll(x),
    ki_pll = PSY.get_ki_pll(x),
)
get_params_metadata(::PSY.KauraPLL) = (
    ω_lp = ParamsMetadata(DEVICE_PARAM, false, false),
    kp_pll = ParamsMetadata(DEVICE_PARAM, false, true),
    ki_pll = ParamsMetadata(DEVICE_PARAM, false, true),
)
#= get_params(x::PSY.FixedFrequency) = [PSY.get_frequency(x)]
get_params_metadata(::PSY.FixedFrequency) =
    [ParamsMetadata(:frequency_FrequencyEstimator, false, false, false, false)]
 =#
#CONVERTER 
get_params(::PSY.AverageConverter) = (;)
get_params_metadata(::PSY.AverageConverter) = (;)

#= get_params(x::PSY.RenewableEnergyConverterTypeA) = [
    PSY.get_T_g(x),
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
    PSY.get_X_source(x),
]
get_params_metadata(::PSY.RenewableEnergyConverterTypeA) = [
    ParamsMetadata(:T_g_Converter, false, false, false, false)
    ParamsMetadata(:Rrpwr_Converter, false, false, false, false)
    ParamsMetadata(:Brkpt_Converter, false, false, false, false)
    ParamsMetadata(:Zerox_Converter, false, false, false, false)
    ParamsMetadata(:Lvpl1_Converter, false, false, false, false)
    ParamsMetadata(:Vo_lim_Converter, false, false, true, false)
    ParamsMetadata(:Lv_pnt0_Converter, false, false, false, false)
    ParamsMetadata(:Lv_pnt1_Converter, false, false, true, false)
    ParamsMetadata(:Io_lim_Converter, false, false, true, false)
    ParamsMetadata(:T_fltr_Converter, false, false, false, false)
    ParamsMetadata(:K_hv_Converter, false, false, false, false)
    ParamsMetadata(:Iqr_min_Converter, false, false, false, false)
    ParamsMetadata(:Iqr_max_Converter, false, false, false, false)
    ParamsMetadata(:Accel_Converter, false, false, false, false)
    ParamsMetadata(:Q_ref_Converter, false, false, false, false)
    ParamsMetadata(:R_source_Converter, false, false, false, false)
    ParamsMetadata(:X_source_Converter, false, false, false, false)
]
 =#
########### GENERATORS #############
function get_params(g::PSY.DynamicGenerator)
    return (
        Machine = get_params(PSY.get_machine(g)),
        Shaft = get_params(PSY.get_shaft(g)),
        AVR = get_params(PSY.get_avr(g)),
        TurbineGov = get_params(PSY.get_prime_mover(g)),
        PSS = get_params(PSY.get_pss(g)),
    )
end
function get_params_metadata(g::PSY.DynamicGenerator)
    return (
        Machine = get_params_metadata(PSY.get_machine(g)),
        Shaft = get_params_metadata(PSY.get_shaft(g)),
        AVR = get_params_metadata(PSY.get_avr(g)),
        TurbineGov = get_params_metadata(PSY.get_prime_mover(g)),
        PSS = get_params_metadata(PSY.get_pss(g)),
    )
end

#MACHINES 
get_params(x::PSY.BaseMachine) =
    (R = PSY.get_R(x), Xd_p = PSY.get_Xd_p(x), eq_p = PSY.get_eq_p(x))
get_params_metadata(::PSY.BaseMachine) = (
    R = ParamsMetadata(DEVICE_PARAM, false, true),
    Xd_p = ParamsMetadata(DEVICE_PARAM, false, true),
    eq_p = ParamsMetadata(DEVICE_PARAM, false, true),
)
get_params(x::PSY.OneDOneQMachine) = (
    R = PSY.get_R(x),
    Xd = PSY.get_Xd(x),
    Xq = PSY.get_Xq(x),
    Xd_p = PSY.get_Xd_p(x),
    Xq_p = PSY.get_Xq_p(x),
    Td0_p = PSY.get_Td0_p(x),
    Tq0_p = PSY.get_Tq0_p(x),
)
get_params_metadata(::PSY.OneDOneQMachine) = (
    R = PartamsMetadata(DEVICE_PARAM, false, true),
    Xd = PartamsMetadata(DEVICE_PARAM, false, true),
    Xq = PartamsMetadata(DEVICE_PARAM, false, true),
    Xd_p = PartamsMetadata(DEVICE_PARAM, false, true),
    Xq_p = PartamsMetadata(DEVICE_PARAM, false, true),
    Td0_p = PartamsMetadata(DEVICE_PARAM, false, false),
    Tq0_p = PartamsMetadata(DEVICE_PARAM, false, false),
)
#= 
#TODO - SimpleMarconatoMachine
get_params(x::PSY.MarconatoMachine) = [
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
    PSY.get_Tq0_pp(x),
    PSY.get_T_AA(x),
    PSY.get_γd(x),
    PSY.get_γq(x),
]
get_params_metadata(::PSY.MarconatoMachine) = [
    ParamsMetadata(:R_Machine, false, false, true, false),
    ParamsMetadata(:Xd_Machine, false, false, true, false),
    ParamsMetadata(:Xq_Machine, false, false, true, false),
    ParamsMetadata(:Xd_p_Machine, false, false, true, false),
    ParamsMetadata(:Xq_p_Machine, false, false, true, false),
    ParamsMetadata(:Xd_pp_Machine, false, false, true, false),
    ParamsMetadata(:Xq_pp_Machine, false, false, true, false),
    ParamsMetadata(:Td0_p_Machine, false, false, true, false),
    ParamsMetadata(:Tq0_p_Machine, false, false, false, false),
    ParamsMetadata(:Td0_pp_Machine, false, false, false, false),
    ParamsMetadata(:Tq0_pp_Machine, false, false, false, false),
    ParamsMetadata(:T_AA_Machine, false, false, true, false),
    ParamsMetadata(:γd_Machine, false, false, true, false),
    ParamsMetadata(:γq_Machine, false, false, true, false),
]
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
get_params_metadata(::PSY.AndersonFouadMachine) = [
    ParamsMetadata(:R_Machine, false, false, true, false),
    ParamsMetadata(:Xd_Machine, false, false, true, false),
    ParamsMetadata(:Xq_Machine, false, false, true, false),
    ParamsMetadata(:Xd_p_Machine, false, false, true, false),
    ParamsMetadata(:Xq_p_Machine, false, false, true, false),
    ParamsMetadata(:Xd_pp_Machine, false, false, true, false),
    ParamsMetadata(:Xq_pp_Machine, false, false, true, false),
    ParamsMetadata(:Td0_p_Machine, false, false, false, false),
    ParamsMetadata(:Tq0_p_Machine, false, false, false, false),
    ParamsMetadata(:Td0_pp_Machine, false, false, false, false),
    ParamsMetadata(:Tq0_pp_Machine, false, false, false, false),
]
#NOTE: Saturation not considered as paramters
=#
get_params(
    x::Union{PSY.RoundRotorMachine, PSY.RoundRotorExponential, PSY.RoundRotorQuadratic},
) = (
    R = PSY.get_R(x),
    Td0_p = PSY.get_Td0_p(x),
    Td0_pp = PSY.get_Td0_pp(x),
    Tq0_p = PSY.get_Tq0_p(x),
    Tq0_pp = PSY.get_Tq0_pp(x),
    Xd = PSY.get_Xd(x),
    Xq = PSY.get_Xq(x),
    Xd_p = PSY.get_Xd_p(x),
    Xq_p = PSY.get_Xq_p(x),
    Xd_pp = PSY.get_Xd_pp(x),
    Xl = PSY.get_Xl(x),
    γ_d1 = PSY.get_γ_d1(x),
    γ_q1 = PSY.get_γ_q1(x),
    γ_d2 = PSY.get_γ_d2(x),
    γ_q2 = PSY.get_γ_q2(x),
    γ_qd = PSY.get_γ_qd(x),
)
get_params_metadata(
    ::Union{PSY.RoundRotorMachine, PSY.RoundRotorExponential, PSY.RoundRotorQuadratic},
) = (
    R = ParamsMetadata(DEVICE_PARAM, false, true),
    Td0_p = ParamsMetadata(DEVICE_PARAM, false, true),
    Td0_pp = ParamsMetadata(DEVICE_PARAM, false, true),
    Tq0_p = ParamsMetadata(DEVICE_PARAM, false, true),
    Tq0_pp = ParamsMetadata(DEVICE_PARAM, false, true),
    Xd = ParamsMetadata(DEVICE_PARAM, false, true),
    Xq = ParamsMetadata(DEVICE_PARAM, false, true),
    Xd_p = ParamsMetadata(DEVICE_PARAM, false, true),
    Xq_p = ParamsMetadata(DEVICE_PARAM, false, true),
    Xd_pp = ParamsMetadata(DEVICE_PARAM, false, true),
    Xl = ParamsMetadata(DEVICE_PARAM, false, true),
    γ_d1 = ParamsMetadata(DEVICE_PARAM, false, true),
    γ_q1 = ParamsMetadata(DEVICE_PARAM, false, true),
    γ_d2 = ParamsMetadata(DEVICE_PARAM, false, true),
    γ_q2 = ParamsMetadata(DEVICE_PARAM, false, true),
    γ_qd = ParamsMetadata(DEVICE_PARAM, false, true),
)
#= get_params(
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
get_params_metadata(
    ::Union{PSY.SalientPoleMachine, PSY.SalientPoleExponential, PSY.SalientPoleQuadratic},
) = [
    ParamsMetadata(:R_Machine, false, false, true, false),
    ParamsMetadata(:Td0_p_Machine, false, false, true, false),
    ParamsMetadata(:Td0_pp_Machine, false, false, true, false),
    ParamsMetadata(:Tq0_pp_Machine, false, false, true, false),
    ParamsMetadata(:Xd_Machine, false, false, true, false),
    ParamsMetadata(:Xq_Machine, false, false, true, false),
    ParamsMetadata(:Xd_p_Machine, false, false, true, false),
    ParamsMetadata(:Xd_pp_Machine, false, false, true, false),
    ParamsMetadata(:Xl_Machine, false, false, true, false),
    ParamsMetadata(:γ_d1_Machine, false, false, true, false),
    ParamsMetadata(:γ_q1_Machine, false, false, true, false),
    ParamsMetadata(:γ_d2_Machine, false, false, true, false),
]  =#

#SHAFTS
get_params(x::PSY.SingleMass) = (H = PSY.get_H(x), D = PSY.get_D(x))
get_params_metadata(::PSY.SingleMass) = (
    H = ParamsMetadata(DEVICE_PARAM, false, false),
    D = ParamsMetadata(DEVICE_PARAM, false, false),
)

#= get_params(x::PSY.FiveMassShaft) = [
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
get_params_metadata(::PSY.FiveMassShaft) = [
    ParamsMetadata(:H_Shaft, false, false, false, false)
    ParamsMetadata(:H_hp_Shaft, false, false, false, false)
    ParamsMetadata(:H_ip_Shaft, false, false, false, false)
    ParamsMetadata(:H_lp_Shaft, false, false, false, false)
    ParamsMetadata(:H_ex_Shaft, false, false, false, false)
    ParamsMetadata(:D_Shaft, false, false, true, false)
    ParamsMetadata(:D_hp_Shaft, false, false, true, false)
    ParamsMetadata(:D_ip_Shaft, false, false, true, false)
    ParamsMetadata(:D_lp_Shaft, false, false, true, false)
    ParamsMetadata(:D_ex_Shaft, false, false, true, false)
    ParamsMetadata(:D_12_Shaft, false, false, true, false)
    ParamsMetadata(:D_23_Shaft, false, false, true, false)
    ParamsMetadata(:D_34_Shaft, false, false, true, false)
    ParamsMetadata(:D_45_Shaft, false, false, true, false)
    ParamsMetadata(:K_hp_Shaft, false, false, true, false)
    ParamsMetadata(:K_ip_Shaft, false, false, true, false)
    ParamsMetadata(:K_lp_Shaft, false, false, true, false)
    ParamsMetadata(:K_ex_Shaft, false, false, true, false)
] =#

#AVRS 
get_params(::PSY.AVRFixed) = (;)
get_params_metadata(::PSY.AVRFixed) = (;)

#= get_params(x::PSY.AVRSimple) = [PSY.get_Kv(x)]
get_params_metadata(::PSY.AVRSimple) = [ParamsMetadata(:Kv_AVR, false, false, false, false)]
 =#
get_params(x::PSY.AVRTypeI) = (
    Ka = PSY.get_Ka(x),
    Ke = PSY.get_Ke(x),
    Kf = PSY.get_Kf(x),
    Ta = PSY.get_Ta(x),
    Te = PSY.get_Te(x),
    Tf = PSY.get_Tf(x),
    Tr = PSY.get_Tr(x),
    Ae = PSY.get_Ae(x),
    Be = PSY.get_Be(x),
)
get_params_metadata(::PSY.AVRTypeI) = (
    Ka = ParamsMetadata(DEVICE_PARAM, false, true),
    Ke = ParamsMetadata(DEVICE_PARAM, false, true),
    Kf = ParamsMetadata(DEVICE_PARAM, false, true),
    Ta = ParamsMetadata(DEVICE_PARAM, false, false),
    Te = ParamsMetadata(DEVICE_PARAM, false, false),
    Tf = ParamsMetadata(DEVICE_PARAM, false, true),
    Tr = ParamsMetadata(DEVICE_PARAM, false, false),
    Ae = ParamsMetadata(DEVICE_PARAM, false, true),
    Be = ParamsMetadata(DEVICE_PARAM, false, true),
)
get_params(x::PSY.SEXS) = (
    Ta_Tb = PSY.get_Ta_Tb(x),
    Tb = PSY.get_Tb(x),
    K = PSY.get_K(x),
    Te = PSY.get_Te(x),
    V_lim = PSY.get_V_lim(x),
)
get_params_metadata(::PSY.SEXS) = (
    Ta_Tb = ParamsMetadata(DEVICE_PARAM, false, true),
    Tb = ParamsMetadata(DEVICE_PARAM, false, false),
    K = ParamsMetadata(DEVICE_PARAM, false, true),
    Te = ParamsMetadata(DEVICE_PARAM, false, false),
    V_lim = (
        min = ParamsMetadata(DEVICE_PARAM, false, true),
        max = ParamsMetadata(DEVICE_PARAM, false, true),
    ),
)
#= get_params(x::PSY.AVRTypeII) = [
    PSY.get_K0(x),
    PSY.get_T1(x),
    PSY.get_T2(x),
    PSY.get_T3(x),
    PSY.get_T4(x),
    PSY.get_Te(x),
    PSY.get_Tr(x),
    PSY.get_Va_lim(x)[1],
    PSY.get_Va_lim(x)[2],
    PSY.get_Ae(x),
    PSY.get_Be(x),
]
get_params_metadata(::PSY.AVRTypeII) = [
    ParamsMetadata(:K0_AVR, false, false, true, false),
    ParamsMetadata(:T1_AVR, false, false, true, false),
    ParamsMetadata(:T2_AVR, false, false, true, false),
    ParamsMetadata(:T3_AVR, false, false, true, false),
    ParamsMetadata(:T4_AVR, false, false, true, false),
    ParamsMetadata(:Te_AVR, false, false, true, false),
    ParamsMetadata(:Tr_AVR, false, false, false, false),
    ParamsMetadata(:Va_min_AVR, false, false, true, false),
    ParamsMetadata(:Va_max_AVR, false, false, true, false),
    ParamsMetadata(:Ae_AVR, false, false, true, false),
    ParamsMetadata(:Be_AVR, false, false, true, false),
]
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
    PSY.get_Vr_lim(x)[2],
]
get_params_metadata(::PSY.ESAC1A) = [
    ParamsMetadata(:Tr, false, false, true, false),
    ParamsMetadata(:Tb_AVR, false, false, true, false),
    ParamsMetadata(:Tc_AVR, false, false, true, false),
    ParamsMetadata(:Ka_AVR, false, false, true, false),
    ParamsMetadata(:Ta_AVR, false, false, false, false),
    ParamsMetadata(:Va_min_AVR, false, false, false, false),
    ParamsMetadata(:Va_max_AVR, false, false, false, false),
    ParamsMetadata(:Te_AVR, false, false, false, false),
    ParamsMetadata(:Kf_AVR, false, false, true, false),
    ParamsMetadata(:Tf_AVR, false, false, true, false),
    ParamsMetadata(:Kc_AVR, false, false, true, false),
    ParamsMetadata(:Kd_AVR, false, false, true, false),
    ParamsMetadata(:Ke_AVR, false, false, true, false),
    ParamsMetadata(:Vr_min_AVR, false, false, true, false),
    ParamsMetadata(:Vr_max_AVR, false, false, true, false),
]
get_params(x::PSY.EXST1) = [
    PSY.get_Tr(x),
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
    PSY.get_Tf(x),
]
get_params_metadata(::PSY.EXST1) = [
    ParamsMetadata(:Tr_AVR, true, false, false, false),
    ParamsMetadata(:Vi_min_AVR, false, false, false, false),
    ParamsMetadata(:Vi_max_AVR, false, false, false, false),
    ParamsMetadata(:Tc_AVR, false, false, true, false),
    ParamsMetadata(:Tb_AVR, true, false, true, false),
    ParamsMetadata(:Ka_AVR, false, false, true, false),
    ParamsMetadata(:Ta_AVR, true, false, false, false),
    ParamsMetadata(:Vr_min_AVR, false, false, true, false),
    ParamsMetadata(:Vr_max_AVR, false, false, true, false),
    ParamsMetadata(:Kc_AVR, false, false, true, false),
    ParamsMetadata(:Kf_AVR, false, false, true, false),
    ParamsMetadata(:Tf_AVR, false, false, true, false),
]
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
get_params_metadata(::PSY.EXAC1) = [
    ParamsMetadata(:Tr_AVR, true, false, false, false),
    ParamsMetadata(:Tb_AVR, true, false, true, false),
    ParamsMetadata(:Tc_AVR, false, false, true, false),
    ParamsMetadata(:Ka_AVR, false, false, true, false),
    ParamsMetadata(:Ta_AVR, true, false, false, false),
    ParamsMetadata(:Vr_min_AVR, false, false, true, false),
    ParamsMetadata(:Vr_max_AVR, false, false, true, false),
    ParamsMetadata(:Te_AVR, false, false, false, false),
    ParamsMetadata(:Kf_AVR, false, false, true, false),
    ParamsMetadata(:Tf_AVR, false, false, true, false),
    ParamsMetadata(:Kc_AVR, false, false, true, false),
    ParamsMetadata(:Kd_AVR, false, false, true, false),
    ParamsMetadata(:Ke_AVR, false, false, true, false),
]
 =#
#TurbineGov
get_params(x::PSY.TGFixed) = (; efficiency = PSY.get_efficiency(x))
get_params_metadata(::PSY.TGFixed) =
    (; efficiency = ParamsMetadata(DEVICE_PARAM, false, true))
get_params(x::PSY.DEGOV) = (
    T1 = PSY.get_T1(x),
    T2 = PSY.get_T2(x),
    T3 = PSY.get_T3(x),
    K = PSY.get_K(x),
    T4 = PSY.get_T4(x),
    T5 = PSY.get_T5(x),
    T6 = PSY.get_T6(x),
    Td = PSY.get_Td(x),
)
get_params_metadata(::PSY.DEGOV) = (
    T1 = ParamsMetadata(DEVICE_PARAM, true, false),
    T2 = ParamsMetadata(DEVICE_PARAM, true, false),
    T3 = ParamsMetadata(DEVICE_PARAM, false, false),
    K = ParamsMetadata(DEVICE_PARAM, false, false),
    T4 = ParamsMetadata(DEVICE_PARAM, false, false),
    T5 = ParamsMetadata(DEVICE_PARAM, true, false),
    T6 = ParamsMetadata(DEVICE_PARAM, true, false),
    Td = ParamsMetadata(DEVICE_PARAM, false, false),
)
#= 
get_params(x::PSY.TGTypeII) = [PSY.get_R(x), PSY.get_T1(x), PSY.get_T2(x)]
get_params_metadata(::PSY.TGTypeII) = [
    ParamsMetadata(:R_tg_TurbineGovR, false, false, true, false),
    ParamsMetadata(:T1_TurbineGov, false, false, true, false),
    ParamsMetadata(:T2_TurbineGov, false, false, true, false),
]
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
get_params_metadata(::PSY.GasTG) = [
    ParamsMetadata(:R_TurbineGov, false, false, true, false),
    ParamsMetadata(:T1_TurbineGov, false, false, false, false),
    ParamsMetadata(:T2_TurbineGov, false, false, false, false),
    ParamsMetadata(:T3_TurbineGov, false, false, false, false),
    ParamsMetadata(:AT_TurbineGov, false, false, true, false),
    ParamsMetadata(:Kt_TurbineGov, false, false, true, false),
    ParamsMetadata(:V_min_TurbineGov, false, false, true, false),
    ParamsMetadata(:V_max_TurbineGov, false, false, true, false),
    ParamsMetadata(:D_turb_TurbineGov, false, false, true, false),
]
get_params(x::PSY.TGTypeI) = [
    PSY.get_R(x),
    PSY.get_Ts(x),
    PSY.get_Tc(x),
    PSY.get_T3(x),
    PSY.get_T4(x),
    PSY.get_T5(x),
    PSY.get_valve_position_limits(x)[1],
    PSY.get_valve_position_limits(x)[2],
]
get_params_metadata(::PSY.TGTypeI) = [
    ParamsMetadata(:R_tg_TurbineGov, false, false, true, false),
    ParamsMetadata(:Ts_TurbineGov, false, false, false, false),
    ParamsMetadata(:Tc_TurbineGov, false, false, true, false),
    ParamsMetadata(:T3_TurbineGov, false, false, true, false),
    ParamsMetadata(:T4_TurbineGov, false, false, true, false),
    ParamsMetadata(:T5_TurbineGov, false, false, true, false),
    ParamsMetadata(:valve_position_min_TurbineGov, false, false, true, false),
    ParamsMetadata(:valve_position_max_TurbineGov, false, false, true, false),
]
 =#

get_params(x::PSY.SteamTurbineGov1) = (
    R = PSY.get_R(x),
    T1 = PSY.get_T1(x),
    valve_position_limits = PSY.get_valve_position_limits(x),
    T2 = PSY.get_T2(x),
    T3 = PSY.get_T3(x),
    D_T = PSY.get_D_T(x),
)
get_params_metadata(::PSY.SteamTurbineGov1) = (
    R = ParamsMetadata(DEVICE_PARAM, false, true),
    T1 = ParamsMetadata(DEVICE_PARAM, false, true),
    valve_position = (
        min = ParamsMetadata(DEVICE_PARAM, false, true),
        max = ParamsMetadata(DEVICE_PARAM, false, true),
    ),
    T2 = ParamsMetadata(DEVICE_PARAM, false, true),
    T3 = ParamsMetadata(DEVICE_PARAM, false, true),
    D_T = ParamsMetadata(DEVICE_PARAM, false, true),
)
#= get_params(x::PSY.HydroTurbineGov) = [
    PSY.get_R(x),
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
    PSY.get_q_nl(x),
]
get_params_metadata(::PSY.HydroTurbineGov) = [
    ParamsMetadata(:R_TurbineGov, false, false, true, false),
    ParamsMetadata(:r_TurbineGov, false, false, true, false),
    ParamsMetadata(:Tr_TurbineGov, false, false, true, false),
    ParamsMetadata(:Tf_TurbineGov, true, false, false, false),
    ParamsMetadata(:Tg_TurbineGov, true, false, false, false),
    ParamsMetadata(:VELM_TurbineGov, false, false, false, false),
    ParamsMetadata(:G_min_TurbineGov, false, false, true, false),
    ParamsMetadata(:G_max_TurbineGov, false, false, true, false),
    ParamsMetadata(:Tw_TurbineGov, false, false, false, false),
    ParamsMetadata(:At_TurbineGov, false, false, true, false),
    ParamsMetadata(:D_T_TurbineGov, false, false, true, false),
    ParamsMetadata(:q_nl_TurbineGov, false, false, true, false),
] =#

#PSS
get_params(x::PSY.PSSFixed) = (; V_pss = PSY.get_V_pss(x))
get_params_metadata(::PSY.PSSFixed) = (; V_pss = ParamsMetadata(DEVICE_PARAM, false, false))
#= 
get_params(x::PSY.STAB1) = [
    PSY.get_KT(x),
    PSY.get_T(x),
    PSY.get_T1T3(x),
    PSY.get_T3(x),
    PSY.get_T2T4(x),
    PSY.get_T4(x),
    PSY.get_H_lim(x),
]
get_params_metadata(::PSY.STAB1) = [
    ParamsMetadata(:KT_PSS, false, false, false, false),
    ParamsMetadata(:T_PSS, false, false, false, false),
    ParamsMetadata(:T1T3_PSS, false, false, false, false),
    ParamsMetadata(:T3_PSS, false, false, false, false),
    ParamsMetadata(:T2T4_PSS, false, false, false, false),
    ParamsMetadata(:T4_PSS, false, false, false, false),
    ParamsMetadata(:H_lim_PSS, false, false, false, false),
] =#

#SOURCE 
get_params(x::PSY.Source) = (
    R_th = PSY.get_R_th(x),
    X_th = PSY.get_X_th(x),
)
get_params_metadata(::PSY.Source) = (
    R_th = ParamsMetadata(DEVICE_PARAM, false, true),
    X_th = ParamsMetadata(DEVICE_PARAM, false, true),
)
#Parameters not implemented for PeriodicVariableSource - requires change in PSY Struct to have information required to construct and deconstruct parameter vector
#= 
#DYNAMIC LOADS 
get_params(x::PSY.ActiveConstantPowerLoad) = [
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
    PSY.get_base_power(x),
]
get_params_metadata(::PSY.ActiveConstantPowerLoad) = [
    ParamsMetadata(:r_load, false, false, true, false),
    ParamsMetadata(:c_dc, false, false, false, false),
    ParamsMetadata(:rf, false, false, true, false),
    ParamsMetadata(:lf, true, false, true, false),
    ParamsMetadata(:cf, true, false, true, false),
    ParamsMetadata(:rg, false, false, true, false),
    ParamsMetadata(:lg, true, false, true, false),
    ParamsMetadata(:kp_pll, false, false, false, false),
    ParamsMetadata(:ki_pll, false, false, false, false),
    ParamsMetadata(:kpv, false, false, false, false),
    ParamsMetadata(:kiv, false, false, true, false),
    ParamsMetadata(:kpc, false, false, false, false),
    ParamsMetadata(:kic, false, false, true, false),
    ParamsMetadata(:base_power, false, false, true, false),
]
get_params(x::PSY.SingleCageInductionMachine) = [
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
get_params_metadata(::PSY.SingleCageInductionMachine) = [
    ParamsMetadata(:R_s, false, false, true, false),
    ParamsMetadata(:R_r, false, false, true, false),
    ParamsMetadata(:X_ls, false, false, true, false),
    ParamsMetadata(:X_lr, false, false, true, false),
    ParamsMetadata(:X_m, false, false, false, false),
    ParamsMetadata(:H, false, false, false, false),
    ParamsMetadata(:A, false, false, true, false),
    ParamsMetadata(:B, false, false, true, false),
    ParamsMetadata(:base_power, false, false, true, false),
    ParamsMetadata(:C, false, false, true, false),
    ParamsMetadata(:τ_ref, false, false, true, false),
    ParamsMetadata(:B_shunt, false, false, true, false),
    ParamsMetadata(:X_ad, false, false, true, false),
    ParamsMetadata(:X_aq, false, false, true, false),
]
get_params(x::PSY.SimplifiedSingleCageInductionMachine) = [
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
get_params_metadata(::PSY.SimplifiedSingleCageInductionMachine) = [
    ParamsMetadata(:R_s, false, false, true, false),
    ParamsMetadata(:R_r, false, false, true, false),
    ParamsMetadata(:X_ls, false, false, false, false),
    ParamsMetadata(:X_lr, false, false, false, false),
    ParamsMetadata(:X_m, false, false, true, false),
    ParamsMetadata(:H, false, false, false, false),
    ParamsMetadata(:A, false, false, true, false),
    ParamsMetadata(:B, false, false, true, false),
    ParamsMetadata(:base_power, false, false, true, false),
    ParamsMetadata(:C, false, false, true, false),
    ParamsMetadata(:τ_ref, false, false, true, false),
    ParamsMetadata(:B_shunt, false, false, true, false),
    ParamsMetadata(:X_ss, false, false, true, false),
    ParamsMetadata(:X_rr, false, false, true, false),
    ParamsMetadata(:X_p, false, false, true, false),
]
get_params(x::PSY.AggregateDistributedGenerationA) = [
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
    PSY.get_Pfa_ref(x),
]
get_params_metadata(::PSY.AggregateDistributedGenerationA) = [
    ParamsMetadata(:T_rv, true, false, false, false),
    ParamsMetadata(:Trf, true, false, false, false),
    ParamsMetadata(:dbd_pnts, false, false, false, false),
    ParamsMetadata(:dbd_pnts, false, false, false, false),
    ParamsMetadata(:K_qv, false, false, true, false),
    ParamsMetadata(:Tp, true, false, false, false),
    ParamsMetadata(:T_iq, true, false, false, false),
    ParamsMetadata(:D_dn, false, false, false, false),
    ParamsMetadata(:D_up, false, false, false, false),
    ParamsMetadata(:fdbd_pnts, false, false, false, false),
    ParamsMetadata(:fdbd_pnts, false, false, false, false),
    ParamsMetadata(:fe_lim, false, false, false, false),
    ParamsMetadata(:fe_lim, false, false, false, false),
    ParamsMetadata(:P_lim, false, false, false, false),
    ParamsMetadata(:P_lim, false, false, false, false),
    ParamsMetadata(:dP_lim, false, false, false, false),
    ParamsMetadata(:dP_lim, false, false, false, false),
    ParamsMetadata(:Tpord, true, false, false, false),
    ParamsMetadata(:Kpg, false, false, false, false),
    ParamsMetadata(:Kig, false, false, false, false),
    ParamsMetadata(:I_max, false, false, false, false),
    ParamsMetadata(:Tg, false, false, false, false),
    ParamsMetadata(:rrpwr, false, false, false, false),
    ParamsMetadata(:Tv, true, false, false, false),
    ParamsMetadata(:Vpr, false, false, false, false),
    ParamsMetadata(:Iq_lim, false, false, false, false),
    ParamsMetadata(:Iq_lim, false, false, false, false),
    ParamsMetadata(:base_power, false, false, false, false),
    ParamsMetadata(:Pfa_ref, false, false, true, false),
]

get_params(x::PSY.CSVGN1) = [
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
    PSY.get_X_th(x),
]
get_params_metadata(::PSY.CSVGN1) = [
    ParamsMetadata(:K, false, false, false, false),
    ParamsMetadata(:T1, false, false, false, false),
    ParamsMetadata(:T2, false, false, false, false),
    ParamsMetadata(:T3, false, false, false, false),
    ParamsMetadata(:T4, false, false, false, false),
    ParamsMetadata(:T5, false, false, false, false),
    ParamsMetadata(:Rmin, false, false, false, false),
    ParamsMetadata(:Vmax, false, false, false, false),
    ParamsMetadata(:Vmin, false, false, false, false),
    ParamsMetadata(:CBase, false, false, false, false),
    ParamsMetadata(:base_power, false, false, false, false),
    ParamsMetadata(:R_th, false, false, false, false),
    ParamsMetadata(:X_th, false, false, false, false),
]
 =#
