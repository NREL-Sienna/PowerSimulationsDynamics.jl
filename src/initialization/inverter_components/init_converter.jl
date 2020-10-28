function initialize_converter!(
    device_states,
    static::PSY.StaticInjection,
    dyn_data::PSY.DynamicInverter{PSY.AverageConverter, O, IC, DC, P, F},
) where {
    O <: PSY.OuterControl,
    IC <: PSY.InnerControl,
    DC <: PSY.DCSource,
    P <: PSY.FrequencyEstimator,
    F <: PSY.Filter,
}

    #Obtain inner variables for component
    md = get_inner_vars(dyn_data)[md_var]
    mq = get_inner_vars(dyn_data)[mq_var]
    Vdc = get_inner_vars(dyn_data)[Vdc_var]
    θ_oc = get_inner_vars(dyn_data)[θ_oc_var]

    #Transform reference frame to grid reference frame
    m_ri = dq_ri(θ_oc + pi / 2) * [md; mq]

    #Update inner_vars
    get_inner_vars(dyn_data)[Vr_cnv_var] = m_ri[R] * Vdc
    get_inner_vars(dyn_data)[Vi_cnv_var] = m_ri[I] * Vdc
end
