function mdl_DCside_ode!( device::PSY.DynamicInverter{C,O,VC,PSY.FixedDCSource,P,F}) where {C <: PSY.Converter,
                                                   O <: PSY.OuterControl,
                                                   VC <: PSY.VSControl,
                                                   P <: PSY.FrequencyEstimator,
                                                   F <: PSY.Filter}

    #Obtain external states inputs for component

    #Obtain inner variables for component

    #Get parameters
    Vdc = device.dc_source.voltage

    #Obtain indices for component w/r to device

    #Define internal states for converter

    #Compute states... none?

    #Update inner_vars
    device.inner_vars[Vdc_var] = Vdc
end
