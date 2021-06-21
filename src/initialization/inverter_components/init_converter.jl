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

    #Obtain inner variables for component
    md = get_inner_vars(dynamic_device)[md_var]
    mq = get_inner_vars(dynamic_device)[mq_var]
    Vdc = get_inner_vars(dynamic_device)[Vdc_var]
    θ_oc = get_inner_vars(dynamic_device)[θ_oc_var]

    #Transform reference frame to grid reference frame
    m_ri = dq_ri(θ_oc + pi / 2) * [md; mq]

    #Update inner_vars
    set_inner_vars!(dynamic_device, Vr_cnv_var, m_ri[R] * Vdc)
    set_inner_vars!(dynamic_device, Vi_cnv_var, m_ri[I] * Vdc)
end
