function initialize_converter!(
    device_states,
    static::PSY.StaticInjection,
    dynamic_device::PSY.DynamicInverter{PSY.AverageConverter, O, IC, DC, P, F},
) where {
    O <: PSY.OuterControl,
    IC <: PSY.InnerControl,
    DC <: PSY.DCSource,
    P <: PSY.FrequencyEstimator,
    F <: PSY.Filter,
}

    #No Update since Vr_cnv_var was updated in Filter using PF solution
    #md = get_inner_vars(dynamic_device)[md_var]
    #mq = get_inner_vars(dynamic_device)[mq_var]
    #Vdc = get_inner_vars(dynamic_device)[Vdc_var]
    #θ_oc = get_inner_vars(dynamic_device)[θ_oc_var]

    #Transform reference frame to grid reference frame
    #m_ri = dq_ri(θ_oc + pi / 2) * [md; mq]

    #Update inner_vars
    #get_inner_vars(dynamic_device)[Vr_cnv_var] = m_ri[R] * Vdc
    #get_inner_vars(dynamic_device)[Vi_cnv_var] = m_ri[I] * Vdc
end

function initialize_converter!(
    device_states,
    static::PSY.StaticInjection,
    dynamic_device::PSY.DynamicInverter{PSY.REGCA1, O, IC, DC, P, F},
) where {
    O <: PSY.OuterControl,
    IC <: PSY.InnerControl,
    DC <: PSY.DCSource,
    P <: PSY.FrequencyEstimator,
    F <: PSY.Filter,
}
    #Get inner vars
    V_t = sqrt(
        get_inner_vars(dynamic_device)[VR_inv_var]^2 +
        get_inner_vars(dynamic_device)[VI_inv_var]^2,
    )
    Iq = get_inner_vars(dynamic_device)[Iq_cnv_var]
    Ip = get_inner_vars(dynamic_device)[Id_cnv_var]

    converter = PSY.get_converter(dynamic_device)
    Io_lim = PSY.get_Io_lim(converter)
    Vo_lim = PSY.get_Vo_lim(converter)

    # Lv_pnt0 is unused in the initialization
    _, Lv_pnt1 = PSY.get_Lv_pnts(converter)

    if (Iq < Io_lim) || (V_t > Vo_lim) || (V_t < Lv_pnt1)
        error("Power flow solution outside of inverter limits. Update parameters.")
    end

    #Update converter states
    converter_ix = get_local_state_ix(dynamic_device, PSY.REGCA1)
    converter_states = @view device_states[converter_ix]
    converter_states[1] = Ip #Ip_cnv
    converter_states[2] = -Iq #Iq_cnv
    converter_states[3] = V_t #Vmeas

    #Update entry to converter currents (inner currents)
    get_inner_vars(dynamic_device)[Id_ic_var] = Ip
    # Multiplied by -1 because PSS/e doesn't do a conjugate.
    get_inner_vars(dynamic_device)[Iq_ic_var] = -Iq
end
