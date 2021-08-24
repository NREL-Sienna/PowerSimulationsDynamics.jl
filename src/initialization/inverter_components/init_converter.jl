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
}

    #Obtain inner variables for component
    md = inner_vars[md_var]
    mq = inner_vars[mq_var]
    Vdc = inner_vars[Vdc_var]
    θ_oc = inner_vars[θ_oc_var]

    #Transform reference frame to grid reference frame
    m_ri = dq_ri(θ_oc + pi / 2) * [md; mq]

    #Update inner_vars
    inner_vars[Vr_cnv_var] = m_ri[R] * Vdc
    inner_vars[Vi_cnv_var] = m_ri[I] * Vdc
end
