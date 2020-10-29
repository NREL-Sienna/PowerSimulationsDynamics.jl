function initialize_DCside!(
    device_states,
    static::PSY.StaticInjection,
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
