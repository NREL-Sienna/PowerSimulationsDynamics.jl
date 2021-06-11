function mass_matrix_DCside_entries!(
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
    @debug "Using default mass matrix entries $DC $(PSY.get_name(dynamic_device))"
end

function mdl_DCside_ode!(
    device_states,
    output_ode,
    f0,
    ω_sys,
    dynamic_device::PSY.DynamicInverter{C, O, IC, PSY.FixedDCSource, P, F},
) where {
    C <: PSY.Converter,
    O <: PSY.OuterControl,
    IC <: PSY.InnerControl,
    P <: PSY.FrequencyEstimator,
    F <: PSY.Filter,
}

    #Update inner_vars
    get_inner_vars(dynamic_device)[Vdc_var] =
        PSY.get_voltage(PSY.get_dc_source(dynamic_device))
end
