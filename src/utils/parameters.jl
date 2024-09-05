
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
function get_params_metadata(x)
    @error "Parameters metadata not yet defined for component: $(typeof(x))"
    (;)
end
#Overloads for refs 
PSY.get_V_ref(x) = 0.0
PSY.get_ω_ref(x) = 0.0
PSY.get_P_ref(x) = 0.0

get_params_metadata(::T) where {T <: PSY.DynamicComponent} = (;)

function get_params(x::PSY.ActivePowerControl)
    @error "Parameters not yet defined for component: $(typeof(x))"
    (;)
end
function get_params_metadata(x::PSY.ActivePowerControl)
    @error "Parameters metadata not yet defined for component: $(typeof(x))"
    (;)
end
function get_params(x::PSY.ReactivePowerControl)
    @error "Parameters not yet defined for component: $(typeof(x))"
    (;)
end
function get_params_metadata(x::PSY.ReactivePowerControl)
    @error "Parameters metadata not yet defined for component: $(typeof(x))"
    (;)
end

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
get_params(x::PSY.ActivePowerPI) = (
    Kp_p = PSY.get_Kp_p(x),
    Ki_p = PSY.get_Ki_p(x),
    ωz = PSY.get_ωz(x),
)
get_params_metadata(::PSY.ActivePowerPI) = (
    Kp_p = ParamsMetadata(DEVICE_PARAM, false, false),
    Ki_p = ParamsMetadata(DEVICE_PARAM, false, true),
    ωz = ParamsMetadata(DEVICE_PARAM, false, false),
)
get_params(x::PSY.ActivePowerDroop) = (Rp = PSY.get_Rp(x), ωz = PSY.get_ωz(x))
get_params_metadata(::PSY.ActivePowerDroop) = (
    Rp = ParamsMetadata(DEVICE_PARAM, false, false),
    ωz = ParamsMetadata(DEVICE_PARAM, false, false),
)
get_params(x::PSY.VirtualInertia) =
    (Ta = PSY.get_Ta(x), kd = PSY.get_kd(x), kω = PSY.get_kω(x))
get_params_metadata(::PSY.VirtualInertia) = (
    Ta = ParamsMetadata(DEVICE_PARAM, false, false),
    kd = ParamsMetadata(DEVICE_PARAM, false, false),
    kω = ParamsMetadata(DEVICE_PARAM, false, false),
)
get_params(x::PSY.ActiveVirtualOscillator) = (k1 = PSY.get_k1(x), ψ = PSY.get_ψ(x))
get_params_metadata(::PSY.ActiveVirtualOscillator) = (
    k1 = ParamsMetadata(DEVICE_PARAM, false, false),
    ψ = ParamsMetadata(DEVICE_PARAM, false, false),
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
get_params(x::PSY.ReactivePowerPI) = (
    Kp_q = PSY.get_Kp_q(x),
    Ki_q = PSY.get_Ki_q(x),
    ωf = PSY.get_ωf(x),
)
get_params_metadata(::PSY.ReactivePowerPI) = (
    Kp_q = ParamsMetadata(DEVICE_PARAM, false, false),
    Ki_q = ParamsMetadata(DEVICE_PARAM, false, true),
    ωf = ParamsMetadata(DEVICE_PARAM, false, false),
)
get_params(x::PSY.ReactivePowerDroop) = (kq = PSY.get_kq(x), ωf = PSY.get_ωf(x))
get_params_metadata(::PSY.ReactivePowerDroop) = (
    kq = ParamsMetadata(DEVICE_PARAM, false, false),
    ωf = ParamsMetadata(DEVICE_PARAM, false, false),
)
get_params(x::PSY.ReactiveVirtualOscillator) = (k2 = PSY.get_k2(x),)
get_params_metadata(::PSY.ReactiveVirtualOscillator) =
    (k2 = ParamsMetadata(DEVICE_PARAM, false, false),)
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
    dbd_pnts1 = PSY.get_dbd_pnts(x)[1],
    dbd_pnts2 = PSY.get_dbd_pnts(x)[2],
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
        dbd_pnts1 = ParamsMetadata(DEVICE_PARAM, false, false),
        dbd_pnts2 = ParamsMetadata(DEVICE_PARAM, false, false),
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
get_params(x::PSY.CurrentModeControl) = (
    kpc = PSY.get_kpc(x),
    kic = PSY.get_kic(x),
    kffv = PSY.get_kffv(x),
)
get_params_metadata(::PSY.CurrentModeControl) = (
    kpc = ParamsMetadata(DEVICE_PARAM, false, true),
    kic = ParamsMetadata(DEVICE_PARAM, false, true),
    kffv = ParamsMetadata(DEVICE_PARAM, false, true),
)
get_params(x::PSY.RECurrentControlB) = (
    Vdip_lim = PSY.get_Vdip_lim(x),
    T_rv = PSY.get_T_rv(x),
    dbd_pnts1 = PSY.get_dbd_pnts(x)[1],
    dbd_pnts2 = PSY.get_dbd_pnts(x)[2],
    K_qv = PSY.get_K_qv(x),
    Iqinj_lim = PSY.get_Iqinj_lim(x),
    V_ref0 = PSY.get_V_ref0(x),
    K_vp = PSY.get_K_vp(x),
    K_vi = PSY.get_K_vi(x),
    T_iq = PSY.get_T_iq(x),
    I_max = PSY.get_I_max(x),
)
get_params_metadata(::PSY.RECurrentControlB) = (
    Vdip_lim = (
        min = ParamsMetadata(DEVICE_PARAM, false, false),
        max = ParamsMetadata(DEVICE_PARAM, false, false),
    ),
    T_rv = ParamsMetadata(DEVICE_PARAM, true, false),
    dbd_pnts1 = ParamsMetadata(DEVICE_PARAM, false, false),
    dbd_pnts2 = ParamsMetadata(DEVICE_PARAM, false, false),
    K_qv = ParamsMetadata(DEVICE_PARAM, false, false),
    Iqinj_lim = (
        min = ParamsMetadata(DEVICE_PARAM, false, false),
        max = ParamsMetadata(DEVICE_PARAM, false, false),
    ),
    V_ref0 = ParamsMetadata(DEVICE_PARAM, false, true),
    K_vp = ParamsMetadata(DEVICE_PARAM, false, false),
    K_vi = ParamsMetadata(DEVICE_PARAM, false, true),
    T_iq = ParamsMetadata(DEVICE_PARAM, true, false),
    I_max = ParamsMetadata(DEVICE_PARAM, false, false),
)

#DC SOURCE 
get_params(x::PSY.FixedDCSource) = (voltage = PSY.get_voltage(x),)
get_params_metadata(::PSY.FixedDCSource) =
    (voltage = ParamsMetadata(DEVICE_PARAM, false, false),)

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
get_params(x::PSY.ReducedOrderPLL) = (
    ω_lp = PSY.get_ω_lp(x),
    kp_pll = PSY.get_kp_pll(x),
    ki_pll = PSY.get_ki_pll(x),
)
get_params_metadata(::PSY.ReducedOrderPLL) = (
    ω_lp = ParamsMetadata(DEVICE_PARAM, false, false),
    kp_pll = ParamsMetadata(DEVICE_PARAM, false, true),
    ki_pll = ParamsMetadata(DEVICE_PARAM, false, true),
)
get_params(x::PSY.FixedFrequency) = (; frequency = PSY.get_frequency(x))
get_params_metadata(::PSY.FixedFrequency) =
    (; frequency = ParamsMetadata(DEVICE_PARAM, false, true))

#CONVERTER 
get_params(::PSY.AverageConverter) = (;)
get_params_metadata(::PSY.AverageConverter) = (;)
get_params(x::PSY.RenewableEnergyConverterTypeA) = (
    T_g = PSY.get_T_g(x),
    Rrpwr = PSY.get_Rrpwr(x),
    Brkpt = PSY.get_Brkpt(x),
    Zerox = PSY.get_Zerox(x),
    Lvpl1 = PSY.get_Lvpl1(x),
    Vo_lim = PSY.get_Vo_lim(x),
    Lv_pnts = PSY.get_Lv_pnts(x),
    Io_lim = PSY.get_Io_lim(x),
    T_fltr = PSY.get_T_fltr(x),
    K_hv = PSY.get_K_hv(x),
    Iqr_lims = PSY.get_Iqr_lims(x),
    Accel = PSY.get_Accel(x),
    Q_ref = PSY.get_Q_ref(x),
    R_source = PSY.get_R_source(x),
    X_source = PSY.get_X_source(x),
)
get_params_metadata(::PSY.RenewableEnergyConverterTypeA) = (
    T_g = ParamsMetadata(DEVICE_PARAM, false, false),
    Rrpwr = ParamsMetadata(DEVICE_PARAM, false, false),
    Brkpt = ParamsMetadata(DEVICE_PARAM, false, false),
    Zerox = ParamsMetadata(DEVICE_PARAM, false, false),
    Lvpl1 = ParamsMetadata(DEVICE_PARAM, false, false),
    Vo_lim = ParamsMetadata(DEVICE_PARAM, false, true),
    Lv_pnts = (
        min = ParamsMetadata(DEVICE_PARAM, false, false),   #confirm
        max = ParamsMetadata(DEVICE_PARAM, false, true), #confirm 
    ),
    Io_lim = ParamsMetadata(DEVICE_PARAM, false, true),
    T_fltr = ParamsMetadata(DEVICE_PARAM, false, false),
    K_hv = ParamsMetadata(DEVICE_PARAM, false, false),
    Iqr_lims = (
        min = ParamsMetadata(DEVICE_PARAM, false, false),
        max = ParamsMetadata(DEVICE_PARAM, false, false),
    ),
    Accel = ParamsMetadata(DEVICE_PARAM, false, false),
    Q_ref = ParamsMetadata(DEVICE_PARAM, false, false),
    R_source = ParamsMetadata(DEVICE_PARAM, false, false),
    X_source = ParamsMetadata(DEVICE_PARAM, false, false),
)
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
    R = ParamsMetadata(DEVICE_PARAM, false, true),
    Xd = ParamsMetadata(DEVICE_PARAM, false, true),
    Xq = ParamsMetadata(DEVICE_PARAM, false, true),
    Xd_p = ParamsMetadata(DEVICE_PARAM, false, true),
    Xq_p = ParamsMetadata(DEVICE_PARAM, false, true),
    Td0_p = ParamsMetadata(DEVICE_PARAM, false, false),
    Tq0_p = ParamsMetadata(DEVICE_PARAM, false, false),
)
get_params(x::PSY.SauerPaiMachine) = (
    R = PSY.get_R(x),
    Xd = PSY.get_Xd(x),
    Xq = PSY.get_Xq(x),
    Xd_p = PSY.get_Xd_p(x),
    Xq_p = PSY.get_Xq_p(x),
    Xd_pp = PSY.get_Xd_pp(x),
    Xq_pp = PSY.get_Xq_pp(x),
    Xl = PSY.get_Xl(x),
    Td0_p = PSY.get_Td0_p(x),
    Tq0_p = PSY.get_Tq0_p(x),
    Td0_pp = PSY.get_Td0_pp(x),
    Tq0_pp = PSY.get_Tq0_pp(x),
    γ_d1 = PSY.get_γ_d1(x),
    γ_q1 = PSY.get_γ_q1(x),
    γ_d2 = PSY.get_γ_d2(x),
    γ_q2 = PSY.get_γ_q2(x),
)
get_params_metadata(::PSY.SauerPaiMachine) = (
    R = ParamsMetadata(DEVICE_PARAM, false, true),
    Xd = ParamsMetadata(DEVICE_PARAM, false, true),
    Xq = ParamsMetadata(DEVICE_PARAM, false, true),
    Xd_p = ParamsMetadata(DEVICE_PARAM, false, true),
    Xq_p = ParamsMetadata(DEVICE_PARAM, false, true),
    Xd_pp = ParamsMetadata(DEVICE_PARAM, false, true),
    Xq_pp = ParamsMetadata(DEVICE_PARAM, false, true),
    Xl = ParamsMetadata(DEVICE_PARAM, false, true),
    Td0_p = ParamsMetadata(DEVICE_PARAM, false, false),
    Tq0_p = ParamsMetadata(DEVICE_PARAM, false, false),
    Td0_pp = ParamsMetadata(DEVICE_PARAM, false, false),
    Tq0_pp = ParamsMetadata(DEVICE_PARAM, false, false),
    γ_d1 = ParamsMetadata(DEVICE_PARAM, false, true),
    γ_q1 = ParamsMetadata(DEVICE_PARAM, false, true),
    γ_d2 = ParamsMetadata(DEVICE_PARAM, false, true),
    γ_q2 = ParamsMetadata(DEVICE_PARAM, false, true),
)
#TODO - SimpleMarconatoMachine
get_params(x::PSY.MarconatoMachine) = (
    R = PSY.get_R(x),
    Xd = PSY.get_Xd(x),
    Xq = PSY.get_Xq(x),
    Xd_p = PSY.get_Xd_p(x),
    Xq_p = PSY.get_Xq_p(x),
    Xd_pp = PSY.get_Xd_pp(x),
    Xq_pp = PSY.get_Xq_pp(x),
    Td0_p = PSY.get_Td0_p(x),
    Tq0_p = PSY.get_Tq0_p(x),
    Td0_pp = PSY.get_Td0_pp(x),
    Tq0_pp = PSY.get_Tq0_pp(x),
    T_AA = PSY.get_T_AA(x),
    γd = PSY.get_γd(x),
    γq = PSY.get_γq(x),
)
get_params_metadata(::PSY.MarconatoMachine) = (
    R = ParamsMetadata(DEVICE_PARAM, false, true),
    Xd = ParamsMetadata(DEVICE_PARAM, false, true),
    Xq = ParamsMetadata(DEVICE_PARAM, false, true),
    Xd_p = ParamsMetadata(DEVICE_PARAM, false, true),
    Xq_p = ParamsMetadata(DEVICE_PARAM, false, true),
    Xd_pp = ParamsMetadata(DEVICE_PARAM, false, true),
    Xq_pp = ParamsMetadata(DEVICE_PARAM, false, true),
    Td0_p = ParamsMetadata(DEVICE_PARAM, false, true),
    Tq0_p = ParamsMetadata(DEVICE_PARAM, false, false),
    Td0_pp = ParamsMetadata(DEVICE_PARAM, false, false),
    Tq0_pp = ParamsMetadata(DEVICE_PARAM, false, false),
    γ_d = ParamsMetadata(DEVICE_PARAM, false, true),
    γ_q = ParamsMetadata(DEVICE_PARAM, false, true),
)
get_params(x::Union{PSY.AndersonFouadMachine, PSY.SimpleAFMachine}) = (
    R = PSY.get_R(x),
    Xd = PSY.get_Xd(x),
    Xq = PSY.get_Xq(x),
    Xd_p = PSY.get_Xd_p(x),
    Xq_p = PSY.get_Xq_p(x),
    Xd_pp = PSY.get_Xd_pp(x),
    Xq_pp = PSY.get_Xq_pp(x),
    Td0_p = PSY.get_Td0_p(x),
    Tq0_p = PSY.get_Tq0_p(x),
    Td0_pp = PSY.get_Td0_pp(x),
    Tq0_pp = PSY.get_Tq0_pp(x),
)
get_params_metadata(::Union{PSY.AndersonFouadMachine, PSY.SimpleAFMachine}) = (
    R = ParamsMetadata(DEVICE_PARAM, false, true),
    Xd = ParamsMetadata(DEVICE_PARAM, false, true),
    Xq = ParamsMetadata(DEVICE_PARAM, false, true),
    Xd_p = ParamsMetadata(DEVICE_PARAM, false, true),
    Xq_p = ParamsMetadata(DEVICE_PARAM, false, true),
    Xd_pp = ParamsMetadata(DEVICE_PARAM, false, true),
    Xq_pp = ParamsMetadata(DEVICE_PARAM, false, true),
    Td0_p = ParamsMetadata(DEVICE_PARAM, false, false),
    Tq0_p = ParamsMetadata(DEVICE_PARAM, false, false),
    Td0_pp = ParamsMetadata(DEVICE_PARAM, false, false),
    Tq0_pp = ParamsMetadata(DEVICE_PARAM, false, false),
)
#NOTE: Saturation not considered as paramters
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

get_params(
    x::Union{PSY.SalientPoleMachine, PSY.SalientPoleExponential, PSY.SalientPoleQuadratic},
) = (
    R = PSY.get_R(x),
    Td0_p = PSY.get_Td0_p(x),
    Td0_pp = PSY.get_Td0_pp(x),
    Tq0_pp = PSY.get_Tq0_pp(x),
    Xd = PSY.get_Xd(x),
    Xq = PSY.get_Xq(x),
    Xd_p = PSY.get_Xd_p(x),
    Xd_pp = PSY.get_Xd_pp(x),
    Xl = PSY.get_Xl(x),
    γ_d1 = PSY.get_γ_d1(x),
    γ_q1 = PSY.get_γ_q1(x),
    γ_d2 = PSY.get_γ_d2(x),
)
get_params_metadata(
    ::Union{PSY.SalientPoleMachine, PSY.SalientPoleExponential, PSY.SalientPoleQuadratic},
) = (
    R = ParamsMetadata(DEVICE_PARAM, false, true),
    Td0_p = ParamsMetadata(DEVICE_PARAM, false, true),
    Td0_pp = ParamsMetadata(DEVICE_PARAM, false, true),
    Tq0_pp = ParamsMetadata(DEVICE_PARAM, false, true),
    Xd = ParamsMetadata(DEVICE_PARAM, false, true),
    Xq = ParamsMetadata(DEVICE_PARAM, false, true),
    Xd_p = ParamsMetadata(DEVICE_PARAM, false, true),
    Xd_pp = ParamsMetadata(DEVICE_PARAM, false, true),
    Xl = ParamsMetadata(DEVICE_PARAM, false, true),
    γ_d1 = ParamsMetadata(DEVICE_PARAM, false, true),
    γ_q1 = ParamsMetadata(DEVICE_PARAM, false, true),
    γ_d2 = ParamsMetadata(DEVICE_PARAM, false, true),
)

#SHAFTS
get_params(x::PSY.SingleMass) = (H = PSY.get_H(x), D = PSY.get_D(x))
get_params_metadata(::PSY.SingleMass) = (
    H = ParamsMetadata(DEVICE_PARAM, false, false),
    D = ParamsMetadata(DEVICE_PARAM, false, false),
)
get_params(x::PSY.FiveMassShaft) = (
    H = PSY.get_H(x),
    H_hp = PSY.get_H_hp(x),
    H_ip = PSY.get_H_ip(x),
    H_lp = PSY.get_H_lp(x),
    H_ex = PSY.get_H_ex(x),
    D = PSY.get_D(x),
    D_hp = PSY.get_D_hp(x),
    D_ip = PSY.get_D_ip(x),
    D_lp = PSY.get_D_lp(x),
    D_ex = PSY.get_D_ex(x),
    D_12 = PSY.get_D_12(x),
    D_23 = PSY.get_D_23(x),
    D_34 = PSY.get_D_34(x),
    D_45 = PSY.get_D_45(x),
    K_hp = PSY.get_K_hp(x),
    K_ip = PSY.get_K_ip(x),
    K_lp = PSY.get_K_lp(x),
    K_ex = PSY.get_K_ex(x),
)
get_params_metadata(::PSY.FiveMassShaft) = (
    H = ParamsMetadata(DEVICE_PARAM, false, false),
    H_hp = ParamsMetadata(DEVICE_PARAM, false, false),
    H_ip = ParamsMetadata(DEVICE_PARAM, false, false),
    H_lp = ParamsMetadata(DEVICE_PARAM, false, false),
    H_ex = ParamsMetadata(DEVICE_PARAM, false, false),
    D = ParamsMetadata(DEVICE_PARAM, false, true),
    D_hp = ParamsMetadata(DEVICE_PARAM, false, true),
    D_ip = ParamsMetadata(DEVICE_PARAM, false, true),
    D_lp = ParamsMetadata(DEVICE_PARAM, false, true),
    D_ex = ParamsMetadata(DEVICE_PARAM, false, true),
    D_12 = ParamsMetadata(DEVICE_PARAM, false, true),
    D_23 = ParamsMetadata(DEVICE_PARAM, false, true),
    D_34 = ParamsMetadata(DEVICE_PARAM, false, true),
    D_45 = ParamsMetadata(DEVICE_PARAM, false, true),
    K_hp = ParamsMetadata(DEVICE_PARAM, false, true),
    K_ip = ParamsMetadata(DEVICE_PARAM, false, true),
    K_lp = ParamsMetadata(DEVICE_PARAM, false, true),
    K_ex = ParamsMetadata(DEVICE_PARAM, false, true),
)

#AVRS 
get_params(::PSY.AVRFixed) = (;)
get_params_metadata(::PSY.AVRFixed) = (;)
get_params(x::PSY.AVRSimple) = (Kv = PSY.get_Kv(x),)
get_params_metadata(::PSY.AVRSimple) = (Kv = ParamsMetadata(DEVICE_PARAM, false, false),)
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
get_params(x::PSY.AVRTypeII) = (
    K0 = PSY.get_K0(x),
    T1 = PSY.get_T1(x),
    T2 = PSY.get_T2(x),
    T3 = PSY.get_T3(x),
    T4 = PSY.get_T4(x),
    Te = PSY.get_Te(x),
    Tr = PSY.get_Tr(x),
    Va_lim = PSY.get_Va_lim(x),
    Ae = PSY.get_Ae(x),
    Be = PSY.get_Be(x),
)
get_params_metadata(::PSY.AVRTypeII) = (
    K0 = ParamsMetadata(DEVICE_PARAM, false, true),
    T1 = ParamsMetadata(DEVICE_PARAM, false, true),
    T2 = ParamsMetadata(DEVICE_PARAM, false, true),
    T3 = ParamsMetadata(DEVICE_PARAM, false, true),
    T4 = ParamsMetadata(DEVICE_PARAM, false, true),
    Te = ParamsMetadata(DEVICE_PARAM, false, true),
    Tr = ParamsMetadata(DEVICE_PARAM, false, false),
    Va_lim = (
        min = ParamsMetadata(DEVICE_PARAM, false, true),
        max = ParamsMetadata(DEVICE_PARAM, false, true),
    ),
    Ae = ParamsMetadata(DEVICE_PARAM, false, true),
    Be = ParamsMetadata(DEVICE_PARAM, false, true),
)
get_params(x::PSY.ESAC1A) = (
    Tr = PSY.get_Tr(x),
    Tb = PSY.get_Tb(x),
    Tc = PSY.get_Tc(x),
    Ka = PSY.get_Ka(x),
    Ta = PSY.get_Ta(x),
    Va_lim = PSY.get_Va_lim(x),
    Te = PSY.get_Te(x),
    Kf = PSY.get_Kf(x),
    Tf = PSY.get_Tf(x),
    Kc = PSY.get_Kc(x),
    Kd = PSY.get_Kd(x),
    Ke = PSY.get_Ke(x),
    Vr_lim = PSY.get_Vr_lim(x),
)
get_params_metadata(::PSY.ESAC1A) = (
    Tr = ParamsMetadata(DEVICE_PARAM, true, false),
    Tb = ParamsMetadata(DEVICE_PARAM, true, false),
    Tc = ParamsMetadata(DEVICE_PARAM, true, false),
    Ka = ParamsMetadata(DEVICE_PARAM, true, false),
    Ta = ParamsMetadata(DEVICE_PARAM, false, false),
    Va_lim = (
        min = ParamsMetadata(DEVICE_PARAM, false, false),
        max = ParamsMetadata(DEVICE_PARAM, false, false),
    ),
    Te = ParamsMetadata(DEVICE_PARAM, false, false),
    Kf = ParamsMetadata(DEVICE_PARAM, true, false),
    Tf = ParamsMetadata(DEVICE_PARAM, true, false),
    Kc = ParamsMetadata(DEVICE_PARAM, true, false),
    Kd = ParamsMetadata(DEVICE_PARAM, true, false),
    Ke = ParamsMetadata(DEVICE_PARAM, true, false),
    Vr_lim = (
        min = ParamsMetadata(DEVICE_PARAM, true, false),
        max = ParamsMetadata(DEVICE_PARAM, true, false),
    ),
)
get_params(x::PSY.EXST1) = (
    Tr = PSY.get_Tr(x),
    Vi_lim = PSY.get_Vi_lim(x),
    Tc = PSY.get_Tc(x),
    Tb = PSY.get_Tb(x),
    Ka = PSY.get_Ka(x),
    Ta = PSY.get_Ta(x),
    Vr_lim = PSY.get_Vr_lim(x),
    Kc = PSY.get_Kc(x),
    Kf = PSY.get_Kf(x),
    Tf = PSY.get_Tf(x),
)
get_params_metadata(::PSY.EXST1) = (
    Tr = ParamsMetadata(DEVICE_PARAM, true, false),
    Vi_lim = (
        min = ParamsMetadata(DEVICE_PARAM, false, false),
        max = ParamsMetadata(DEVICE_PARAM, false, false),
    ),
    Tc = ParamsMetadata(DEVICE_PARAM, false, true),
    Tb = ParamsMetadata(DEVICE_PARAM, true, true),
    Ka = ParamsMetadata(DEVICE_PARAM, false, true),
    Ta = ParamsMetadata(DEVICE_PARAM, true, false),
    Vr_lim = (
        min = ParamsMetadata(DEVICE_PARAM, false, true),
        max = ParamsMetadata(DEVICE_PARAM, false, true),
    ),
    Kc = ParamsMetadata(DEVICE_PARAM, false, true),
    Kf = ParamsMetadata(DEVICE_PARAM, false, true),
    Tf = ParamsMetadata(DEVICE_PARAM, false, true),
)
get_params(x::PSY.EXAC1) = (
    Tr = PSY.get_Tr(x),
    Tb = PSY.get_Tb(x),
    Tc = PSY.get_Tc(x),
    Ka = PSY.get_Ka(x),
    Ta = PSY.get_Ta(x),
    Vr_lim = PSY.get_Vr_lim(x),
    Te = PSY.get_Te(x),
    Kf = PSY.get_Kf(x),
    Tf = PSY.get_Tf(x),
    Kc = PSY.get_Kc(x),
    Kd = PSY.get_Kd(x),
    Ke = PSY.get_Ke(x),
)
get_params_metadata(::PSY.EXAC1) = (
    Tr = ParamsMetadata(DEVICE_PARAM, true, false),
    Tb = ParamsMetadata(DEVICE_PARAM, true, true),
    Tc = ParamsMetadata(DEVICE_PARAM, false, true),
    Ka = ParamsMetadata(DEVICE_PARAM, false, true),
    Ta = ParamsMetadata(DEVICE_PARAM, true, false),
    Vr_lim = (
        min = ParamsMetadata(DEVICE_PARAM, false, true),
        max = ParamsMetadata(DEVICE_PARAM, false, true),
    ),
    Te = ParamsMetadata(DEVICE_PARAM, false, false),
    Kf = ParamsMetadata(DEVICE_PARAM, false, true),
    Tf = ParamsMetadata(DEVICE_PARAM, false, true),
    Kc = ParamsMetadata(DEVICE_PARAM, false, true),
    Kd = ParamsMetadata(DEVICE_PARAM, false, true),
    Ke = ParamsMetadata(DEVICE_PARAM, false, true),
)
get_params(x::PSY.SCRX) = (
    Ta_Tb = PSY.get_Ta_Tb(x),
    Tb = PSY.get_Tb(x),
    K = PSY.get_K(x),
    Te = PSY.get_Te(x),
    Efd_lim = PSY.get_Efd_lim(x),
    rc_rfd = PSY.get_rc_rfd(x),
)
get_params_metadata(::PSY.SCRX) = (
    Ta_Tb = ParamsMetadata(DEVICE_PARAM, false, true),
    Tb = ParamsMetadata(DEVICE_PARAM, true, true),
    K = ParamsMetadata(DEVICE_PARAM, false, true),
    Te = ParamsMetadata(DEVICE_PARAM, true, false),
    Efd_lim = (
        min = ParamsMetadata(DEVICE_PARAM, false, true),
        max = ParamsMetadata(DEVICE_PARAM, false, true),
    ),
    rc_rfd = ParamsMetadata(DEVICE_PARAM, false, true),
)
get_params(x::PSY.ESST1A) = (
    Tr = PSY.get_Tr(x),
    Vi_lim = (min = PSY.get_Vi_lim(x)[1], max = min = PSY.get_Vi_lim(x)[2]),
    Tc = PSY.get_Tc(x),
    Tb = PSY.get_Tb(x),
    Tc1 = PSY.get_Tc1(x),
    Tb1 = PSY.get_Tb1(x),
    Ka = PSY.get_Ka(x),
    Ta = PSY.get_Ta(x),
    Va_lim = PSY.get_Va_lim(x),
    Vr_lim = PSY.get_Vr_lim(x),
    Kc = PSY.get_Kc(x),
    Kf = PSY.get_Kf(x),
    Tf = PSY.get_Tf(x),
    K_lr = PSY.get_K_lr(x),
    I_lr = PSY.get_I_lr(x),
)
get_params_metadata(::PSY.ESST1A) = (
    Tr = ParamsMetadata(DEVICE_PARAM, true, false),
    Vi_lim = (
        min = ParamsMetadata(DEVICE_PARAM, false, false),
        max = ParamsMetadata(DEVICE_PARAM, false, false),
    ),
    Tc = ParamsMetadata(DEVICE_PARAM, false, true),
    Tb = ParamsMetadata(DEVICE_PARAM, true, true),
    Tc1 = ParamsMetadata(DEVICE_PARAM, false, true),
    Tb1 = ParamsMetadata(DEVICE_PARAM, true, true),
    Ka = ParamsMetadata(DEVICE_PARAM, false, true),
    Ta = ParamsMetadata(DEVICE_PARAM, true, false),
    Va_lim = (
        min = ParamsMetadata(DEVICE_PARAM, false, false),
        max = ParamsMetadata(DEVICE_PARAM, false, false),
    ),
    Vr_lim = (
        min = ParamsMetadata(DEVICE_PARAM, false, true),
        max = ParamsMetadata(DEVICE_PARAM, false, true),
    ),
    Kc = ParamsMetadata(DEVICE_PARAM, false, true),
    Kf = ParamsMetadata(DEVICE_PARAM, false, true),
    Tf = ParamsMetadata(DEVICE_PARAM, false, true),
    K_lr = ParamsMetadata(DEVICE_PARAM, false, true),
    I_lr = ParamsMetadata(DEVICE_PARAM, false, true),
)
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
)
get_params_metadata(::PSY.DEGOV) = (
    T1 = ParamsMetadata(DEVICE_PARAM, true, false),
    T2 = ParamsMetadata(DEVICE_PARAM, true, false),
    T3 = ParamsMetadata(DEVICE_PARAM, false, false),
    K = ParamsMetadata(DEVICE_PARAM, false, false),
    T4 = ParamsMetadata(DEVICE_PARAM, false, false),
    T5 = ParamsMetadata(DEVICE_PARAM, true, false),
    T6 = ParamsMetadata(DEVICE_PARAM, true, false),
)
get_params(x::PSY.TGTypeII) = (
    R = PSY.get_R(x),
    T1 = PSY.get_T1(x),
    T2 = PSY.get_T2(x),
)
get_params_metadata(::PSY.TGTypeII) = (
    R = ParamsMetadata(DEVICE_PARAM, false, true),
    T1 = ParamsMetadata(DEVICE_PARAM, false, true),
    T2 = ParamsMetadata(DEVICE_PARAM, false, true),
)
get_params(x::PSY.GasTG) = (
    R = PSY.get_R(x),
    T1 = PSY.get_T1(x),
    T2 = PSY.get_T2(x),
    T3 = PSY.get_T3(x),
    AT = PSY.get_AT(x),
    Kt = PSY.get_Kt(x),
    V_lim = (min = PSY.get_V_lim(x)[1], max = PSY.get_V_lim(x)[2]),
    D_turb = PSY.get_D_turb(x),
)
get_params_metadata(::PSY.GasTG) = (
    R = ParamsMetadata(DEVICE_PARAM, false, true),
    T1 = ParamsMetadata(DEVICE_PARAM, false, false),
    T2 = ParamsMetadata(DEVICE_PARAM, false, false),
    T3 = ParamsMetadata(DEVICE_PARAM, false, false),
    AT = ParamsMetadata(DEVICE_PARAM, false, true),
    Kt = ParamsMetadata(DEVICE_PARAM, false, true),
    V_lim = (
        min = ParamsMetadata(DEVICE_PARAM, false, true),
        max = ParamsMetadata(DEVICE_PARAM, false, true),
    ),
    D_turb = ParamsMetadata(DEVICE_PARAM, false, true),
)

get_params(x::PSY.TGTypeI) = (
    R = PSY.get_R(x),
    Ts = PSY.get_Ts(x),
    Tc = PSY.get_Tc(x),
    T3 = PSY.get_T3(x),
    T4 = PSY.get_T4(x),
    T5 = PSY.get_T5(x),
    valve_position_limits = PSY.get_valve_position_limits(x),
)

get_params_metadata(::PSY.TGTypeI) = (
    R = ParamsMetadata(DEVICE_PARAM, false, true),
    Ts = ParamsMetadata(DEVICE_PARAM, false, false),
    Tc = ParamsMetadata(DEVICE_PARAM, false, true),
    T3 = ParamsMetadata(DEVICE_PARAM, false, true),
    T4 = ParamsMetadata(DEVICE_PARAM, false, true),
    T5 = ParamsMetadata(DEVICE_PARAM, false, true),
    valve_position_limits = (
        min = ParamsMetadata(DEVICE_PARAM, false, true),
        max = ParamsMetadata(DEVICE_PARAM, false, true),
    ),
)
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
    valve_position_limits = (
        min = ParamsMetadata(DEVICE_PARAM, false, true),
        max = ParamsMetadata(DEVICE_PARAM, false, true),
    ),
    T2 = ParamsMetadata(DEVICE_PARAM, false, true),
    T3 = ParamsMetadata(DEVICE_PARAM, false, true),
    D_T = ParamsMetadata(DEVICE_PARAM, false, true),
)
get_params(x::PSY.HydroTurbineGov) = (
    R = PSY.get_R(x),
    r = PSY.get_r(x),
    Tr = PSY.get_Tr(x),
    Tf = PSY.get_Tf(x),
    Tg = PSY.get_Tg(x),
    VELM = PSY.get_VELM(x),
    gate_position_limits = PSY.get_gate_position_limits(x),
    Tw = PSY.get_Tw(x),
    At = PSY.get_At(x),
    D_T = PSY.get_D_T(x),
    q_nl = PSY.get_q_nl(x),
)
get_params_metadata(::PSY.HydroTurbineGov) = (
    R = ParamsMetadata(DEVICE_PARAM, false, true),
    r = ParamsMetadata(DEVICE_PARAM, false, true),
    Tr = ParamsMetadata(DEVICE_PARAM, false, true),
    Tf = ParamsMetadata(DEVICE_PARAM, true, false),
    Tg = ParamsMetadata(DEVICE_PARAM, true, false),
    VELM = ParamsMetadata(DEVICE_PARAM, false, false),
    gate_position_limits = (
        min = ParamsMetadata(DEVICE_PARAM, false, true),
        max = ParamsMetadata(DEVICE_PARAM, false, true),
    ),
    Tw = ParamsMetadata(DEVICE_PARAM, false, false),
    At = ParamsMetadata(DEVICE_PARAM, false, true),
    D_T = ParamsMetadata(DEVICE_PARAM, false, true),
    q_nl = ParamsMetadata(DEVICE_PARAM, false, true),
)
#PSS
get_params(x::PSY.PSSFixed) = (; V_pss = PSY.get_V_pss(x))
get_params_metadata(::PSY.PSSFixed) = (; V_pss = ParamsMetadata(DEVICE_PARAM, false, false))
get_params(x::PSY.STAB1) = (
    KT = PSY.get_KT(x),
    T = PSY.get_T(x),
    T1T3 = PSY.get_T1T3(x),
    T3 = PSY.get_T3(x),
    T2T4 = PSY.get_T2T4(x),
    T4 = PSY.get_T4(x),
    H_lim = PSY.get_H_lim(x),
)
get_params_metadata(::PSY.STAB1) = (
    KT = ParamsMetadata(DEVICE_PARAM, false, false),
    T = ParamsMetadata(DEVICE_PARAM, false, false),
    T1T3 = ParamsMetadata(DEVICE_PARAM, false, false),
    T3 = ParamsMetadata(DEVICE_PARAM, false, false),
    T2T4 = ParamsMetadata(DEVICE_PARAM, false, false),
    T4 = ParamsMetadata(DEVICE_PARAM, false, false),
    H_lim = ParamsMetadata(DEVICE_PARAM, false, false),
)
get_params(x::PSY.IEEEST) = (
    A1 = PSY.get_A1(x),
    A2 = PSY.get_A2(x),
    A3 = PSY.get_A3(x),
    A4 = PSY.get_A4(x),
    A5 = PSY.get_A5(x),
    A6 = PSY.get_A6(x),
    T1 = PSY.get_T1(x),
    T2 = PSY.get_T2(x),
    T3 = PSY.get_T3(x),
    T4 = PSY.get_T4(x),
    T5 = PSY.get_T5(x),
    T6 = PSY.get_T6(x),
    Ks = PSY.get_Ks(x),
    Ls_lim1 = PSY.get_Ls_lim(x)[1],
    Ls_lim2 = PSY.get_Ls_lim(x)[2],
    Vcu = PSY.get_Vcu(x),
    Vcl = PSY.get_Vcl(x),
)
get_params_metadata(::PSY.IEEEST) = (
    A1 = ParamsMetadata(DEVICE_PARAM, false, true),
    A2 = ParamsMetadata(DEVICE_PARAM, true, true),
    A3 = ParamsMetadata(DEVICE_PARAM, false, false),
    A4 = ParamsMetadata(DEVICE_PARAM, true, false),
    A5 = ParamsMetadata(DEVICE_PARAM, false, true),
    A6 = ParamsMetadata(DEVICE_PARAM, false, true),
    T1 = ParamsMetadata(DEVICE_PARAM, false, true),
    T2 = ParamsMetadata(DEVICE_PARAM, true, true),
    T3 = ParamsMetadata(DEVICE_PARAM, false, true),
    T4 = ParamsMetadata(DEVICE_PARAM, true, true),
    T5 = ParamsMetadata(DEVICE_PARAM, false, true),
    T6 = ParamsMetadata(DEVICE_PARAM, true, true),
    Ks = ParamsMetadata(DEVICE_PARAM, false, true),
    Ls_lim1 = ParamsMetadata(DEVICE_PARAM, false, true),
    Ls_lim2 = ParamsMetadata(DEVICE_PARAM, false, true),
    Vcu = ParamsMetadata(DEVICE_PARAM, false, true),
    Vcl = ParamsMetadata(DEVICE_PARAM, false, true),
)
get_params(x::PSY.PSS2B) = (
    Tw1 = PSY.get_Tw1(x),
    Tw2 = PSY.get_Tw2(x),
    T6 = PSY.get_T6(x),
    Tw3 = PSY.get_Tw3(x),
    Tw4 = PSY.get_Tw4(x),
    T7 = PSY.get_T7(x),
    Ks2 = PSY.get_Ks2(x),
    Ks3 = PSY.get_Ks3(x),
    T8 = PSY.get_T8(x),
    T9 = PSY.get_T9(x),
    Ks1 = PSY.get_Ks1(x),
    T1 = PSY.get_T1(x),
    T2 = PSY.get_T2(x),
    T3 = PSY.get_T3(x),
    T4 = PSY.get_T4(x),
    T10 = PSY.get_T10(x),
    T11 = PSY.get_T11(x),
)
get_params_metadata(::PSY.PSS2B) = (
    Tw1 = ParamsMetadata(DEVICE_PARAM, false, true),
    Tw2 = ParamsMetadata(DEVICE_PARAM, false, false),
    T6 = ParamsMetadata(DEVICE_PARAM, false, false),
    Tw3 = ParamsMetadata(DEVICE_PARAM, false, true),
    Tw4 = ParamsMetadata(DEVICE_PARAM, false, false),
    T7 = ParamsMetadata(DEVICE_PARAM, false, false),
    Ks2 = ParamsMetadata(DEVICE_PARAM, false, false),
    Ks3 = ParamsMetadata(DEVICE_PARAM, false, false),
    T8 = ParamsMetadata(DEVICE_PARAM, false, false),
    T9 = ParamsMetadata(DEVICE_PARAM, false, true),
    Ks1 = ParamsMetadata(DEVICE_PARAM, false, false),
    T1 = ParamsMetadata(DEVICE_PARAM, false, false),
    T2 = ParamsMetadata(DEVICE_PARAM, false, false),
    T3 = ParamsMetadata(DEVICE_PARAM, false, false),
    T4 = ParamsMetadata(DEVICE_PARAM, false, false),
    T10 = ParamsMetadata(DEVICE_PARAM, false, false),
    T11 = ParamsMetadata(DEVICE_PARAM, false, false),
)
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

#DYNAMIC LOADS 
get_params(x::PSY.ActiveConstantPowerLoad) = (
    r_load = PSY.get_r_load(x),
    c_dc = PSY.get_c_dc(x),
    rf = PSY.get_rf(x),
    lf = PSY.get_lf(x),
    cf = PSY.get_cf(x),
    rg = PSY.get_rg(x),
    lg = PSY.get_lg(x),
    kp_pll = PSY.get_kp_pll(x),
    ki_pll = PSY.get_ki_pll(x),
    kpv = PSY.get_kpv(x),
    kiv = PSY.get_kiv(x),
    kpc = PSY.get_kpc(x),
    kic = PSY.get_kic(x),
    base_power = PSY.get_base_power(x),
)
get_params_metadata(::PSY.ActiveConstantPowerLoad) = (
    r_load = ParamsMetadata(DEVICE_PARAM, false, true),
    c_dc = ParamsMetadata(DEVICE_PARAM, false, false),
    rf = ParamsMetadata(DEVICE_PARAM, false, true),
    lf = ParamsMetadata(DEVICE_PARAM, true, true),
    cf = ParamsMetadata(DEVICE_PARAM, true, true),
    rg = ParamsMetadata(DEVICE_PARAM, false, true),
    lg = ParamsMetadata(DEVICE_PARAM, true, true),
    kp_pll = ParamsMetadata(DEVICE_PARAM, false, false),
    ki_pll = ParamsMetadata(DEVICE_PARAM, false, false),
    kpv = ParamsMetadata(DEVICE_PARAM, false, false),
    kiv = ParamsMetadata(DEVICE_PARAM, false, true),
    kpc = ParamsMetadata(DEVICE_PARAM, false, false),
    kic = ParamsMetadata(DEVICE_PARAM, false, true),
    base_power = ParamsMetadata(DEVICE_PARAM, false, false),
)
get_params(x::PSY.SingleCageInductionMachine) = (
    R_s = PSY.get_R_s(x),
    R_r = PSY.get_R_r(x),
    X_ls = PSY.get_X_ls(x),
    X_lr = PSY.get_X_lr(x),
    X_m = PSY.get_X_m(x),
    H = PSY.get_H(x),
    A = PSY.get_A(x),
    B = PSY.get_B(x),
    base_power = PSY.get_base_power(x),
    C = PSY.get_C(x),
    τ_ref = PSY.get_τ_ref(x),
    B_shunt = PSY.get_B_shunt(x),
    X_ad = PSY.get_X_ad(x),
    X_aq = PSY.get_X_aq(x),
)
get_params_metadata(::PSY.SingleCageInductionMachine) = (
    R_s = ParamsMetadata(DEVICE_PARAM, false, true),
    R_r = ParamsMetadata(DEVICE_PARAM, false, true),
    X_ls = ParamsMetadata(DEVICE_PARAM, false, true),
    X_lr = ParamsMetadata(DEVICE_PARAM, false, true),
    X_m = ParamsMetadata(DEVICE_PARAM, false, false),
    H = ParamsMetadata(DEVICE_PARAM, false, false),
    A = ParamsMetadata(DEVICE_PARAM, false, true),
    B = ParamsMetadata(DEVICE_PARAM, false, true),
    base_power = ParamsMetadata(DEVICE_PARAM, false, true),
    C = ParamsMetadata(DEVICE_PARAM, false, true),
    τ_ref = ParamsMetadata(DEVICE_PARAM, false, true),
    B_shunt = ParamsMetadata(DEVICE_PARAM, false, true),
    X_ad = ParamsMetadata(DEVICE_PARAM, false, true),
    X_aq = ParamsMetadata(DEVICE_PARAM, false, true),
)

get_params(x::PSY.SimplifiedSingleCageInductionMachine) = (
    R_s = PSY.get_R_s(x),
    R_r = PSY.get_R_r(x),
    X_ls = PSY.get_X_ls(x),
    X_lr = PSY.get_X_lr(x),
    X_m = PSY.get_X_m(x),
    H = PSY.get_H(x),
    A = PSY.get_A(x),
    B = PSY.get_B(x),
    base_power = PSY.get_base_power(x),
    C = PSY.get_C(x),
    τ_ref = PSY.get_τ_ref(x),
    B_shunt = PSY.get_B_shunt(x),
    X_ss = PSY.get_X_ss(x),
    X_rr = PSY.get_X_rr(x),
    X_p = PSY.get_X_p(x),
)
get_params_metadata(::PSY.SimplifiedSingleCageInductionMachine) = (
    R_s = ParamsMetadata(DEVICE_PARAM, false, true),
    R_r = ParamsMetadata(DEVICE_PARAM, false, true),
    X_ls = ParamsMetadata(DEVICE_PARAM, false, false),
    X_lr = ParamsMetadata(DEVICE_PARAM, false, false),
    X_m = ParamsMetadata(DEVICE_PARAM, false, true),
    H = ParamsMetadata(DEVICE_PARAM, false, false),
    A = ParamsMetadata(DEVICE_PARAM, false, true),
    B = ParamsMetadata(DEVICE_PARAM, false, true),
    base_power = ParamsMetadata(DEVICE_PARAM, false, true),
    C = ParamsMetadata(DEVICE_PARAM, false, true),
    τ_ref = ParamsMetadata(DEVICE_PARAM, false, true),
    B_shunt = ParamsMetadata(DEVICE_PARAM, false, true),
    X_ss = ParamsMetadata(DEVICE_PARAM, false, true),
    X_rr = ParamsMetadata(DEVICE_PARAM, false, true),
    X_p = ParamsMetadata(DEVICE_PARAM, false, true),
)
get_params(x::PSY.AggregateDistributedGenerationA) = (
    T_rv = PSY.get_T_rv(x),
    Trf = PSY.get_Trf(x),
    dbd_pnts1 = PSY.get_dbd_pnts(x)[1],
    dbd_pnts2 = PSY.get_dbd_pnts(x)[2],
    K_qv = PSY.get_K_qv(x),
    Tp = PSY.get_Tp(x),
    T_iq = PSY.get_T_iq(x),
    D_dn = PSY.get_D_dn(x),
    D_up = PSY.get_D_up(x),
    fdbd_pnts1 = PSY.get_fdbd_pnts(x)[1],
    fdbd_pnts2 = PSY.get_fdbd_pnts(x)[2],
    fe_lim = PSY.get_fe_lim(x),
    P_lim = PSY.get_P_lim(x),
    dP_lim = PSY.get_dP_lim(x),
    Tpord = PSY.get_Tpord(x),
    Kpg = PSY.get_Kpg(x),
    Kig = PSY.get_Kig(x),
    I_max = PSY.get_I_max(x),
    Tg = PSY.get_Tg(x),
    rrpwr = PSY.get_rrpwr(x),
    Tv = PSY.get_Tv(x),
    Vpr = PSY.get_Vpr(x),
    Iq_lim = PSY.get_Iq_lim(x),
    base_power = PSY.get_base_power(x),
    Pfa_ref = PSY.get_Pfa_ref(x),
)
get_params_metadata(::PSY.AggregateDistributedGenerationA) = (
    T_rv = ParamsMetadata(DEVICE_PARAM, true, false),
    Trf = ParamsMetadata(DEVICE_PARAM, true, false),
    dbd_pnts1 = ParamsMetadata(DEVICE_PARAM, false, false),
    dbd_pnts2 = ParamsMetadata(DEVICE_PARAM, false, false),
    K_qv = ParamsMetadata(DEVICE_PARAM, false, true),
    Tp = ParamsMetadata(DEVICE_PARAM, true, false),
    T_iq = ParamsMetadata(DEVICE_PARAM, true, false),
    D_dn = ParamsMetadata(DEVICE_PARAM, false, false),
    D_up = ParamsMetadata(DEVICE_PARAM, false, false),
    fdbd_pnts1 = ParamsMetadata(DEVICE_PARAM, false, false),
    fdbd_pnts2 = ParamsMetadata(DEVICE_PARAM, false, false),
    fe_lim = (
        min = ParamsMetadata(DEVICE_PARAM, false, false),
        max = ParamsMetadata(DEVICE_PARAM, false, false),
    ),
    P_lim = (
        min = ParamsMetadata(DEVICE_PARAM, false, false),
        max = ParamsMetadata(DEVICE_PARAM, false, false),
    ),
    dP_lim = (
        min = ParamsMetadata(DEVICE_PARAM, false, false),
        max = ParamsMetadata(DEVICE_PARAM, false, false),
    ),
    Tpord = ParamsMetadata(DEVICE_PARAM, true, false),
    Kpg = ParamsMetadata(DEVICE_PARAM, false, false),
    Kig = ParamsMetadata(DEVICE_PARAM, false, false),
    I_max = ParamsMetadata(DEVICE_PARAM, false, false),
    Tg = ParamsMetadata(DEVICE_PARAM, false, false),
    rrpwr = ParamsMetadata(DEVICE_PARAM, false, false),
    Tv = ParamsMetadata(DEVICE_PARAM, true, false),
    Vpr = ParamsMetadata(DEVICE_PARAM, false, false),
    Iq_lim = (
        min = ParamsMetadata(DEVICE_PARAM, false, false),
        max = ParamsMetadata(DEVICE_PARAM, false, false),
    ),
    base_power = ParamsMetadata(DEVICE_PARAM, false, false),
    Pfa_ref = ParamsMetadata(DEVICE_PARAM, false, true),
)
get_params(x::PSY.CSVGN1) = (
    K = PSY.get_K(x),
    T1 = PSY.get_T1(x),
    T2 = PSY.get_T2(x),
    T3 = PSY.get_T3(x),
    T4 = PSY.get_T4(x),
    T5 = PSY.get_T5(x),
    Rmin = PSY.get_Rmin(x),
    Vmax = PSY.get_Vmax(x),
    Vmin = PSY.get_Vmin(x),
    CBase = PSY.get_CBase(x),
    base_power = PSY.get_base_power(x),
    R_th = PSY.get_R_th(x),
    X_th = PSY.get_X_th(x),
)
get_params_metadata(::PSY.CSVGN1) = (
    K = ParamsMetadata(DEVICE_PARAM, false, false),
    T1 = ParamsMetadata(DEVICE_PARAM, false, false),
    T2 = ParamsMetadata(DEVICE_PARAM, false, false),
    T3 = ParamsMetadata(DEVICE_PARAM, false, false),
    T4 = ParamsMetadata(DEVICE_PARAM, false, false),
    T5 = ParamsMetadata(DEVICE_PARAM, false, false),
    Rmin = ParamsMetadata(DEVICE_PARAM, false, false),
    Vmax = ParamsMetadata(DEVICE_PARAM, false, false),
    Vmin = ParamsMetadata(DEVICE_PARAM, false, false),
    CBase = ParamsMetadata(DEVICE_PARAM, false, false),
    base_power = ParamsMetadata(DEVICE_PARAM, false, false),
    R_th = ParamsMetadata(DEVICE_PARAM, false, false),
    X_th = ParamsMetadata(DEVICE_PARAM, false, false),
)
