function mdl_converter_ode!(
    device::PSY.DynamicInverter{PSY.AvgCnvFixedDC,O,VC,DC,P,F},
) where {
    O<:PSY.OuterControl,
    VC<:PSY.VSControl,
    DC<:PSY.DCSource,
    P<:PSY.FrequencyEstimator,
    F<:PSY.Filter,
}

    #Obtain external states inputs for component

    #Obtain inner variables for component
    md = get_inner_vars(device)[md_var]
    mq = get_inner_vars(device)[mq_var]
    VDC = get_inner_vars(device)[Vdc_var]

    #Update inner_vars
    get_inner_vars(device)[Vdcnv_var] = md * VDC
    get_inner_vars(device)[Vqcnv_var] = mq * VDC
end

#TODO: Same as above, but:
## VDC is state of DC source
## IDC is algebraic inner_var output from Vcnv real output power
