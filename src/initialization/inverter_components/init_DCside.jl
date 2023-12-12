function initialize_DCside!(
    device_states,
    static::PSY.StaticInjection,
    dynamic_device::DynamicWrapper{
        PSY.DynamicInverter{C, O, IC, PSY.FixedDCSource, P, F, L},
    },
    inner_vars::AbstractVector,
) where {
    C <: PSY.Converter,
    O <: PSY.OuterControl,
    IC <: PSY.InnerControl,
    P <: PSY.FrequencyEstimator,
    F <: PSY.Filter,
    L <: Union{Nothing, PSY.InverterLimiter},
}

    #Update inner_vars
    inner_vars[Vdc_var] = PSY.get_voltage(PSY.get_dc_source(dynamic_device))
    return
end
