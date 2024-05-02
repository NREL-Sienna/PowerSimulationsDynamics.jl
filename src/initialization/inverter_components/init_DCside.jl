function initialize_DCside!(
    device_states,
    device_parameters,
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
    local_ix_params = get_local_parameter_ix(dynamic_device, PSY.FixedDCSource)
    internal_params = @view device_parameters[local_ix_params]
    #Update inner_vars
    inner_vars[Vdc_var] = internal_params[1]
    return
end
