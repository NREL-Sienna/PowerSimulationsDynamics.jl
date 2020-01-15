function mdl_converter_ode!( device::PSY.DynamicInverter{PSY.AvgCnvFixedDC,O,VC,DC,P,F}) where {O <: PSY.OuterControl,
                                                   VC<: PSY.VSControl,
                                                   DC<: PSY.DCSource,
                                                   P <: PSY.FrequencyEstimator,
                                                   F <: PSY.Filter}

    #Obtain external states inputs for component

    #Obtain inner variables for component
     md = device.inner_vars[md_var]
     mq = device.inner_vars[mq_var]
    VDC = device.inner_vars[Vdc_var]

    #Get parameters

    #Obtain indices for component w/r to device

    #Define internal states for converter

    #Compute states... none?

    #Update inner_vars
    device.inner_vars[Vdcnv_var] = md*VDC
    device.inner_vars[Vqcnv_var] = mq*VDC
    #device.inner_vars[Vdcnv_var] = 1.02
    #device.inner_vars[Vqcnv_var] = 0.01
end

#TODO: Same as above, but:
## VDC is state of DC source
## IDC is algebraic inner_var output from Vcnv real output power
