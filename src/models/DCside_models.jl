function mdl_DCside_ode!( device::DynInverter{C,O,VC,FixedDCSource,P,F}) where {C <: Converter,
                                                   O <: OuterControl,
                                                   VC <: VSControl,
                                                   P <: FrequencyEstimator,
                                                   F <: Filter}

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
