function mdl_DCside_ode!(
    dyn_data::PSY.DynamicInverter{C, O, IC, PSY.FixedDCSource, P, F},
) where {
    C <: PSY.Converter,
    O <: PSY.OuterControl,
    IC <: PSY.InnerControl,
    P <: PSY.FrequencyEstimator,
    F <: PSY.Filter,
}

    #Update inner_vars
    get_inner_vars(dyn_data)[Vdc_var] = PSY.get_voltage(PSY.get_dc_source(dyn_data))
end
