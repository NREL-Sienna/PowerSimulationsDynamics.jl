function initialize_converter!(
    device_states,
    static::PSY.StaticInjection,
    dynamic_device::PSY.DynamicInverter{PSY.AverageConverter, O, IC, DC, P, F},
    inner_vars::AbstractVector,
) where {
    O <: PSY.OuterControl,
    IC <: PSY.InnerControl,
    DC <: PSY.DCSource,
    P <: PSY.FrequencyEstimator,
    F <: PSY.Filter,
} end

function initialize_converter!(
    device_states,
    static::PSY.StaticInjection,
    dynamic_device::PSY.DynamicInverter{PSY.RenewableEnergyConverterTypeA, O, IC, DC, P, F},
) where {
    O <: PSY.OuterControl,
    IC <: PSY.InnerControl,
    DC <: PSY.DCSource,
    P <: PSY.FrequencyEstimator,
    F <: PSY.Filter,
}
    #Get inner vars
    V_R = get_inner_vars(dynamic_device)[Vr_inv_var]
    V_I = get_inner_vars(dynamic_device)[Vi_inv_var]
    θ = atan(V_I / V_R)
    V_t = sqrt(V_R^2 + V_I^2)
    Iq_external = get_inner_vars(dynamic_device)[Ii_cnv_var]
    Ip_external = get_inner_vars(dynamic_device)[Ir_cnv_var]
    #Reference Transformation
    Ip = Ip_external * cos(-θ) - Iq_external * sin(-θ)
    Iq = Ip_external * sin(-θ) + Iq_external * cos(-θ)
    converter = PSY.get_converter(dynamic_device)
    Io_lim = PSY.get_Io_lim(converter)
    Vo_lim = PSY.get_Vo_lim(converter)

    # Lv_pnt0 is unused in the initialization
    _, Lv_pnt1 = PSY.get_Lv_pnts(converter)

    #Obtain inner variables for component
    md = inner_vars[md_var]
    mq = inner_vars[mq_var]
    Vdc = inner_vars[Vdc_var]
    θ_oc = inner_vars[θ_oc_var]

    if (Iq < Io_lim) || (V_t > Vo_lim) || (V_t < Lv_pnt1)
        error("Power flow solution outside of inverter limits. Update parameters.")
    end

    #Update converter states
    converter_ix = get_local_state_ix(dynamic_device, PSY.RenewableEnergyConverterTypeA)
    converter_states = @view device_states[converter_ix]
    converter_states[1] = Ip #Ip_cnv
    converter_states[2] = -Iq #Iq_cnv
    converter_states[3] = V_t #Vmeas

    #Update inner_vars
    inner_vars[Vr_cnv_var] = m_ri[R] * Vdc
    inner_vars[Vi_cnv_var] = m_ri[I] * Vdc
end
