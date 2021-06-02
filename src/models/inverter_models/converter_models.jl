function mass_matrix_converter_entries!(
    mass_matrix,
    dynamic_device::PSY.DynamicInverter{C, O, IC, DC, P, F},
) where {
    C <: PSY.Converter,
    O <: PSY.OuterControl,
    IC <: PSY.InnerControl,
    DC <: PSY.DCSource,
    P <: PSY.FrequencyEstimator,
    F <: PSY.Filter,
}
    @debug "Using default mass matrix entries $C $(PSY.get_name(dynamic_device))"
end

function mdl_converter_ode!(
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
