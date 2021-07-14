function mass_matrix_converter_entries!(
    mass_matrix,
    converter::C,
    global_index::Dict{Symbol, Int64},
) where {C <: PSY.Converter}
    @debug "Using default mass matrix entries $C"
end

function mdl_converter_ode!(
    device_states,
    output_ode,
    Q_gen0::Float64,
    dynamic_device::PSY.DynamicInverter{PSY.AverageConverter, O, IC, DC, P, F},
) where {
    O <: PSY.OuterControl,
    IC <: PSY.InnerControl,
    DC <: PSY.DCSource,
    P <: PSY.FrequencyEstimator,
    F <: PSY.Filter,
}

    #Obtain external states inputs for component

    #Obtain inner variables for component
    md = get_inner_vars(dynamic_device)[md_var]
    mq = get_inner_vars(dynamic_device)[mq_var]
    Vdc = get_inner_vars(dynamic_device)[Vdc_var]
    θ_oc = get_inner_vars(dynamic_device)[θ_oc_var]

    #Transform reference frame to grid reference frame
    m_ri = dq_ri(θ_oc + pi / 2) * [md; mq]

    #Update inner_vars
    get_inner_vars(dynamic_device)[Vr_cnv_var] = m_ri[R] * Vdc
    get_inner_vars(dynamic_device)[Vi_cnv_var] = m_ri[I] * Vdc
end

function mdl_converter_ode!(
    device_states,
    output_ode,
    dynamic_device::PSY.DynamicInverter{PSY.REGCA1, O, IC, DC, P, F},
) where {
    O <: PSY.OuterControl,
    IC <: PSY.InnerControl,
    DC <: PSY.DCSource,
    P <: PSY.FrequencyEstimator,
    F <: PSY.Filter,
}

    #Obtain external states inputs for component
    #external_ix = get_input_port_ix(dynamic_device, PSY.REGCA1)

    #Obtain inner variables for component
    V_t = sqrt(
        get_inner_vars(dynamic_device)[VR_inv_var]^2 +
        get_inner_vars(dynamic_device)[VI_inv_var]^2,
    )
    Ip_cmd = get_inner_vars(dynamic_device)[Id_ic_var]
    Iq_cmd = get_inner_vars(dynamic_device)[Iq_ic_var]

    #Get Converter parameters
    converter = PSY.get_converter(dynamic_device)
    T_g = PSY.get_T_g(converter)
    Rrpwr = PSY.get_Rrpwr(converter)
    Brkpt = PSY.get_Brkpt(converter)
    Zerox = PSY.get_Zerox(converter)
    Lvpl1 = PSY.get_Lvpl1(converter)
    Vo_lim = PSY.get_Vo_lim(converter)
    Lv_pnt0, Lv_pnt1 = PSY.get_Lv_pnts(converter)
    Io_lim = PSY.get_Io_lim(converter)
    T_fltr = PSY.get_T_fltr(converter)
    K_hv = PSY.get_K_hv(converter)
    Iqr_min, Iqr_max = PSY.get_Iqr_lims(converter)
    #Accel = PSY.get_Accel(converter)
    Lvpl_sw = PSY.get_Lvpl_sw(converter)
    Q_ref = PSY.get_Q_ref(converter)

    #Obtain indices for component w/r to device
    local_ix = get_local_state_ix(dynamic_device, PSY.REGCA1)
    #Define internal states for Converter
    internal_states = @view device_states[local_ix]
    Ip = internal_states[1]
    Iq = internal_states[2]
    Vmeas = internal_states[3]

    # Compute additional variables
    # Active Power Part
    Rp_dn = -Inf
    Rp_up = Inf
    #Update Ip ramp limits
    Ip >= 0.0 ? Rp_up = Rrpwr : Rp_dn = -Rrpwr
    #Saturate Ip if LVPL is active
    Ip_sat = Ip
    Ip_binary = 1.0
    if Lvpl_sw == 1
        LVPL = get_LVPL_gain(Vmeas, Zerox, Brkpt, Lvpl1)
        Ip_sat = Ip <= LVPL ? Ip : LVPL
        Ip_binary = Ip <= LVPL ? 1.0 : 0.0
    end
    Ip_in = clamp(Ip_cmd - Ip_sat, Rp_dn, Rp_up)
    #Get Low Voltage Active Current Management Gain
    G_lv = get_LV_current_gain(V_t, Lv_pnt0, Lv_pnt1)

    # Reactive Power Part
    Rq_dn = -Inf
    Rq_up = Inf
    #Update Iq ramp limits
    Q_ref >= 0 ? Rq_up = Iqr_max : Rq_dn = Iqr_min
    Iq_in = clamp(Iq_cmd - Iq, Rq_dn, Rq_up)
    Iq_extra = max(K_hv * (V_t - Vo_lim), 0.0)

    #Update ODEs
    output_ode[local_ix[1]] = Ip_binary * (1.0 / T_g) * Ip_in #(Ip_cmd - Ip)
    output_ode[local_ix[2]] = (1.0 / T_g) * Iq_in # (Iq_cmd - Iq)
    output_ode[local_ix[3]] = (1.0 / T_fltr) * (V_t - Vmeas)

    #Update inner_vars
    get_inner_vars(dynamic_device)[Id_cnv_var] = G_lv * Ip_sat
    get_inner_vars(dynamic_device)[Iq_cnv_var] = max(-Iq - Iq_extra, Io_lim) #-Iq
end
