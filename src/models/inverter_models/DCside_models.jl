function mdl_DCside_ode!( device::PSY.DynamicInverter{C,O,VC,PSY.FixedDCSource,P,F}) where {C <: PSY.Converter,
                                                   O <: PSY.OuterControl,
                                                   VC <: PSY.VSControl,
                                                   P <: PSY.FrequencyEstimator,
                                                   F <: PSY.Filter}

    #Update inner_vars
    device.inner_vars[Vdc_var] = PSY.get_voltage(PSY.get_dc_source(device))
end
