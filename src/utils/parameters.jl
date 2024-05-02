struct ParamsMetadata
    symbol::Symbol
    in_mass_matrix::Bool
    in_network::Bool
    impacts_ic::Bool
    impacts_pf::Bool
end
get_params_symbol(x) = [metadata.symbol for metadata in get_params_metadata(x)]

# TODO - temporary for dynamic components that have not yet been modified to use parameters.
function get_params(x::PSY.Device)
    @warn "Parameters not yet defined for device: $(typeof(x))"
    Float64[]
end
function get_params(x::T) where {T <: PSY.DynamicComponent}
    @warn "Parameters not yet defined for dynamic component: $(typeof(x))"
    Float64[]
end
get_params_metadata(::T) where {T <: PSY.DynamicComponent} = ParamsMetadata[]
get_params(::PSY.ActivePowerControl) = Float64[]
get_params_metadata(::PSY.ActivePowerControl) = ParamsMetadata[]
get_params(::PSY.ReactivePowerControl) = Float64[]
get_params_metadata(::PSY.ReactivePowerControl) = ParamsMetadata[]

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

get_params(x::PSY.Line) = [PSY.get_r(x), PSY.get_x(x)]
get_params_metadata(::PSY.Line) = [
    ParamsMetadata(:r, false, true, true, true),
    ParamsMetadata(:x, false, true, true, true),
]
get_params(::StaticLoadWrapper) = Float64[]
get_params_metadata(::StaticLoadWrapper) = ParamsMetadata[]
########### INVERTERS #############
function get_params(g::PSY.DynamicInverter)
    vcat(
        get_params(PSY.get_converter(g)),
        get_params(PSY.get_outer_control(g)),
        get_params(PSY.get_inner_control(g)),
        get_params(PSY.get_dc_source(g)),
        get_params(PSY.get_freq_estimator(g)),
        get_params(PSY.get_filter(g)),
    )
end
function get_params_metadata(g::PSY.DynamicInverter)
    vcat(
        get_params_metadata(PSY.get_converter(g)),
        get_params_metadata(PSY.get_outer_control(g)),
        get_params_metadata(PSY.get_inner_control(g)),
        get_params_metadata(PSY.get_dc_source(g)),
        get_params_metadata(PSY.get_freq_estimator(g)),
        get_params_metadata(PSY.get_filter(g)),
    )
end

#FILTERS 
get_params(x::PSY.LCLFilter) =
    [PSY.get_lf(x), PSY.get_rf(x), PSY.get_cf(x), PSY.get_lg(x), PSY.get_rg(x)]
get_params_metadata(::PSY.LCLFilter) = [
    ParamsMetadata(:lf_Filter, true, false, true, false),
    ParamsMetadata(:rf_Filter, false, false, true, false),
    ParamsMetadata(:cf_Filter, true, false, true, false),
    ParamsMetadata(:lg_Filter, true, false, true, false),
    ParamsMetadata(:rg_Filter, false, false, true, false),
]
get_params(x::PSY.RLFilter) = [PSY.get_rf(x), PSY.get_lf(x)]
get_params_metadata(::PSY.RLFilter) = [
    ParamsMetadata(:rf_Filter, false, false, true, false),
    ParamsMetadata(:lf_Filter, false, false, true, false),
]

#OUTER CONTROL
get_params(x::PSY.OuterControl) = vcat(
    get_params(PSY.get_active_power_control(x)),
    get_params(PSY.get_reactive_power_control(x)),
)
get_params_metadata(x::PSY.OuterControl) = vcat(
    get_params_metadata(PSY.get_active_power_control(x)),
    get_params_metadata(PSY.get_reactive_power_control(x)),
)
#ACTIVE POWER CONTROL
get_params(x::PSY.VirtualInertia) = [PSY.get_Ta(x), PSY.get_kd(x), PSY.get_kω(x)]
get_params_metadata(::PSY.VirtualInertia) = [
    ParamsMetadata(:Ta_ActivePowerControl, false, false, false, false),
    ParamsMetadata(:kd_ActivePowerControl, false, false, false, false),
    ParamsMetadata(:kω_ActivePowerControl, false, false, false, false),
]
get_params(x::PSY.ActiveRenewableControllerAB) = [
    PSY.get_K_pg(x),
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
get_params_metadata(::PSY.ActiveRenewableControllerAB) = [
    ParamsMetadata(:K_pg_ActivePowerControl, false, false, false, false),
    ParamsMetadata(:K_ig_ActivePowerControl, false, false, true, false),
    ParamsMetadata(:T_p_ActivePowerControl, false, false, true, false),
    ParamsMetadata(:fdbd1_ActivePowerControl, false, false, false, false),
    ParamsMetadata(:fdbd2_ActivePowerControl, false, false, false, false),
    ParamsMetadata(:fe_min_ActivePowerControl, false, false, false, false),
    ParamsMetadata(:fe_max_ActivePowerControl, false, false, false, false),
    ParamsMetadata(:P_min_ActivePowerControl, false, false, false, false),
    ParamsMetadata(:P_max_ActivePowerControl, false, false, false, false),
    ParamsMetadata(:T_g_ap_ActivePowerControl, false, false, false, false),
    ParamsMetadata(:D_dn_ActivePowerControl, false, false, false, false),
    ParamsMetadata(:D_up_ActivePowerControl, false, false, false, false),
    ParamsMetadata(:dP_min_ActivePowerControl, false, false, false, false),
    ParamsMetadata(:dP_max_ActivePowerControl, false, false, false, false),
    ParamsMetadata(:P_min_inner_ActivePowerControl, false, false, false, false),
    ParamsMetadata(:P_max_inner_ActivePowerControl, false, false, false, false),
    ParamsMetadata(:T_pord_ActivePowerControl, false, false, false, false),
]

#REACTIVE POWER CONTROL
get_params(x::PSY.ReactivePowerDroop) = [PSY.get_kq(x), PSY.get_ωf(x)]
get_params_metadata(x::PSY.ReactivePowerDroop) = [
    ParamsMetadata(:kq_ReactivePowerControl, false, false, false, false),
    ParamsMetadata(:ωf_ReactivePowerControl, false, false, false, false),
]
get_params(x::PSY.ReactiveRenewableControllerAB) = [
    PSY.get_T_fltr(x),
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
    PSY.get_K_qi(x),
]
get_params_metadata(::PSY.ReactiveRenewableControllerAB) = [
    ParamsMetadata(:T_fltr_ReactivePowerControl, false, false, false, false),
    ParamsMetadata(:K_p_ReactivePowerControl, false, false, false, false),
    ParamsMetadata(:K_i_ReactivePowerControl, false, false, true, false),
    ParamsMetadata(:T_ft_ReactivePowerControl, false, false, false, false),
    ParamsMetadata(:T_fv_ReactivePowerControl, false, false, false, false),
    ParamsMetadata(:V_frz_ReactivePowerControl, false, false, false, false),
    ParamsMetadata(:R_c_ReactivePowerControl, false, false, true, false),
    ParamsMetadata(:X_c_ReactivePowerControl, false, false, true, false),
    ParamsMetadata(:K_c_ReactivePowerControl, false, false, true, false),
    ParamsMetadata(:e_min_ReactivePowerControl, false, false, false, false),
    ParamsMetadata(:e_max_ReactivePowerControl, false, false, false, false),
    ParamsMetadata(:dbd_pnts1_ReactivePowerControl, false, false, false, false),
    ParamsMetadata(:dbd_pnts2_ReactivePowerControl, false, false, false, false),
    ParamsMetadata(:Q_min_ReactivePowerControl, false, false, false, false),
    ParamsMetadata(:Q_max_ReactivePowerControl, false, false, false, false),
    ParamsMetadata(:T_p_ReactivePowerControl, false, false, true, false),
    ParamsMetadata(:Q_min_inner_ReactivePowerControl, false, false, false, false),
    ParamsMetadata(:Q_max_inner_ReactivePowerControl, false, false, false, false),
    ParamsMetadata(:V_min_ReactivePowerControl, false, false, false, false),
    ParamsMetadata(:V_max_ReactivePowerControl, false, false, false, false),
    ParamsMetadata(:K_qp_ReactivePowerControl, false, false, false, false),
    ParamsMetadata(:K_qi_ReactivePowerControl, false, false, true, false),
]

#INNER CONTROL
get_params(x::PSY.VoltageModeControl) = [
    PSY.get_kpv(x),
    PSY.get_kiv(x),
    PSY.get_kffv(x),
    PSY.get_rv(x),
    PSY.get_lv(x),
    PSY.get_kpc(x),
    PSY.get_kic(x),
    PSY.get_kffi(x),
    PSY.get_ωad(x),
    PSY.get_kad(x),
]
get_params_metadata(::PSY.VoltageModeControl) = [
    ParamsMetadata(:kpv_InnerControl, false, false, true, false),
    ParamsMetadata(:kiv_InnerControl, false, false, true, false),
    ParamsMetadata(:kffv_InnerControl, false, false, true, false),
    ParamsMetadata(:rv_InnerControl, false, false, true, false),
    ParamsMetadata(:lv_InnerControl, false, false, true, false),
    ParamsMetadata(:kpc_InnerControl, false, false, true, false),
    ParamsMetadata(:kic_InnerControl, false, false, true, false),
    ParamsMetadata(:kffi_InnerControl, false, false, true, false),
    ParamsMetadata(:ωad_InnerControl, false, false, false, false),
    ParamsMetadata(:kad_InnerControl, false, false, true, false),
]
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

#DC SOURCE 
get_params(x::PSY.FixedDCSource) = [PSY.get_voltage(x)]
get_params_metadata(::PSY.FixedDCSource) =
    [ParamsMetadata(:voltage_DCSource, false, false, false, false)]

#FREQ ESTIMATOR
get_params(x::PSY.KauraPLL) = [PSY.get_ω_lp(x), PSY.get_kp_pll(x), PSY.get_ki_pll(x)]
get_params_metadata(::PSY.KauraPLL) = [
    ParamsMetadata(:ω_lp_FrequencyEstimator, false, false, false, false)
    ParamsMetadata(:kp_pll_FrequencyEstimator, false, false, true, false)
    ParamsMetadata(:ki_pll_FrequencyEstimator, false, false, true, false)
]
get_params(x::PSY.FixedFrequency) = [PSY.get_frequency(x)]
get_params_metadata(::PSY.FixedFrequency) =
    [ParamsMetadata(:frequency_FrequencyEstimator, false, false, false, false)]

#CONVERTER 
get_params(::PSY.AverageConverter) = Float64[]
get_params_metadata(::PSY.AverageConverter) = ParamsMetadata[]
get_params(x::PSY.RenewableEnergyConverterTypeA) = [
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

########### GENERATORS #############
function get_params(g::PSY.DynamicGenerator)
    vcat(
        get_params(PSY.get_machine(g)),
        get_params(PSY.get_shaft(g)),
        get_params(PSY.get_avr(g)),
        get_params(PSY.get_prime_mover(g)),
        get_params(PSY.get_pss(g)),
    )
end
function get_params_metadata(g::PSY.DynamicGenerator)
    vcat(
        get_params_metadata(PSY.get_machine(g)),
        get_params_metadata(PSY.get_shaft(g)),
        get_params_metadata(PSY.get_avr(g)),
        get_params_metadata(PSY.get_prime_mover(g)),
        get_params_metadata(PSY.get_pss(g)),
    )
end

#MACHINES 
get_params(x::PSY.BaseMachine) = [PSY.get_R(x), PSY.get_Xd_p(x), PSY.get_eq_p(x)]
get_params_metadata(::PSY.BaseMachine) = [
    ParamsMetadata(:R_Machine, false, false, true, false),
    ParamsMetadata(:Xd_p_Machine, false, false, true, false),
    ParamsMetadata(:eq_p_Machine, false, false, true, false),
]
get_params(x::PSY.OneDOneQMachine) = [
    PSY.get_R(x),
    PSY.get_Xd(x),
    PSY.get_Xq(x),
    PSY.get_Xd_p(x),
    PSY.get_Xq_p(x),
    PSY.get_Td0_p(x),
    PSY.get_Tq0_p(x),
]
get_params_metadata(::PSY.OneDOneQMachine) = [
    ParamsMetadata(:R_Machine, false, false, true, false),
    ParamsMetadata(:Xd_Machine, false, false, true, false),
    ParamsMetadata(:Xq_Machine, false, false, true, false),
    ParamsMetadata(:Xd_p_Machine, false, false, true, false),
    ParamsMetadata(:Xq_p_Machine, false, false, true, false),
    ParamsMetadata(:Td0_p_Machine, false, false, false, false),
    ParamsMetadata(:Tq0_p_Machine, false, false, false, false),
]
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
get_params_metadata(
    ::Union{PSY.RoundRotorMachine, PSY.RoundRotorExponential, PSY.RoundRotorQuadratic},
) = [
    ParamsMetadata(:R_Machine, false, false, true, false),
    ParamsMetadata(:Td0_p_Machine, false, false, true, false),
    ParamsMetadata(:Td0_pp_Machine, false, false, true, false),
    ParamsMetadata(:Tq0_p_Machine, false, false, true, false),
    ParamsMetadata(:Tq0_pp_Machine, false, false, true, false),
    ParamsMetadata(:Xd_Machine, false, false, true, false),
    ParamsMetadata(:Xq_Machine, false, false, true, false),
    ParamsMetadata(:Xd_p_Machine, false, false, true, false),
    ParamsMetadata(:Xq_p_Machine, false, false, true, false),
    ParamsMetadata(:Xd_pp_Machine, false, false, true, false),
    ParamsMetadata(:Xl_Machine, false, false, true, false),
    ParamsMetadata(:γ_d1_Machine, false, false, true, false),
    ParamsMetadata(:γ_q1_Machine, false, false, true, false),
    ParamsMetadata(:γ_d2_Machine, false, false, true, false),
    ParamsMetadata(:γ_q2_Machine, false, false, true, false),
    ParamsMetadata(:γ_qd_Machine, false, false, true, false),
]
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
]

#SHAFTS
get_params(x::PSY.SingleMass) = [PSY.get_H(x), PSY.get_D(x)]
get_params_metadata(::PSY.SingleMass) = [
    ParamsMetadata(:H_Shaft, false, false, false, false),
    ParamsMetadata(:D_Shaft, false, false, false, false),
]
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
]

#AVRS 
get_params(::PSY.AVRFixed) = Float64[]
get_params_metadata(::PSY.AVRFixed) = ParamsMetadata[]
get_params(x::PSY.AVRSimple) = [PSY.get_Kv(x)]
get_params_metadata(::PSY.AVRSimple) = [ParamsMetadata(:Kv_AVR, false, false, false, false)]
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
get_params_metadata(::PSY.AVRTypeI) = [
    ParamsMetadata(:Ka_AVR, false, false, true, false),
    ParamsMetadata(:Ke_AVR, false, false, true, false),
    ParamsMetadata(:Kf_AVR, false, false, true, false),
    ParamsMetadata(:Ta_AVR, false, false, false, false),
    ParamsMetadata(:Te_AVR, false, false, false, false),
    ParamsMetadata(:Tf_AVR, false, false, true, false),
    ParamsMetadata(:Tr_AVR, false, false, false, false),
    ParamsMetadata(:Ae_AVR, false, false, true, false),
    ParamsMetadata(:Be_AVR, false, false, true, false),
]
get_params(x::PSY.SEXS) = [
    PSY.get_Ta_Tb(x),
    PSY.get_Tb(x),
    PSY.get_K(x),
    PSY.get_Te(x),
    PSY.get_V_lim(x)[1],
    PSY.get_V_lim(x)[2],
]
get_params_metadata(::PSY.SEXS) = [
    ParamsMetadata(:Ta_Tb_AVR, false, false, true, false)
    ParamsMetadata(:Tb_AVR, false, false, false, false)
    ParamsMetadata(:K_AVR, false, false, true, false)
    ParamsMetadata(:Te_AVR, false, false, false, false)
    ParamsMetadata(:V_min_AVR, false, false, true, false)
    ParamsMetadata(:V_max_AVR, false, false, true, false)
]
get_params(x::PSY.AVRTypeII) = [
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

#TurbineGov
get_params(x::PSY.TGFixed) = [PSY.get_efficiency(x)]
get_params_metadata(::PSY.TGFixed) =
    [ParamsMetadata(:efficiency_TurbineGov, false, false, true, false)]
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
get_params(x::PSY.SteamTurbineGov1) = [PSY.get_R(x),
    PSY.get_T1(x),
    PSY.get_valve_position_limits(x)[1],
    PSY.get_valve_position_limits(x)[2],
    PSY.get_T2(x),
    PSY.get_T3(x),
    PSY.get_D_T(x)]
get_params_metadata(::PSY.SteamTurbineGov1) = [
    ParamsMetadata(:R_tg, false, false, true, false),
    ParamsMetadata(:T1_TurbineGov, false, false, true, false),
    ParamsMetadata(:valve_position_min_TurbineGov, false, false, true, false),
    ParamsMetadata(:valve_position_max_TurbineGov, false, false, true, false),
    ParamsMetadata(:T2_TurbineGov, false, false, true, false),
    ParamsMetadata(:T3_TurbineGov, false, false, true, false),
    ParamsMetadata(:D_T_TurbineGov, false, false, true, false),
]
get_params(x::PSY.HydroTurbineGov) = [
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
]

#PSS
get_params(x::PSY.PSSFixed) = [PSY.get_V_pss(x)]
get_params_metadata(::PSY.PSSFixed) =
    [ParamsMetadata(:V_pss_PSS, false, false, false, false)]
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
]

#SOURCE 
get_params(x::PSY.Source) = [
    PSY.get_R_th(x),
    PSY.get_X_th(x),
]
get_params_metadata(::PSY.Source) = [
    ParamsMetadata(:R_th, false, false, true, false),
    ParamsMetadata(:X_th, false, false, true, false),
]
#Parameters not implemented for PeriodicVariableSource - requires change in PSY Struct to have information required to construct and deconstruct parameter vector

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
