function mdl_filter_ode!(
    device_states,
    output_ode,
    current_r,
    current_i,
    sys_Sbase,
    f0,
    ω_sys,
    device::PSY.DynamicInverter{C, O, IC, DC, P, PSY.LCLFilter},
) where {
    C <: PSY.Converter,
    O <: PSY.OuterControl,
    IC <: PSY.InnerControl,
    DC <: PSY.DCSource,
    P <: PSY.FrequencyEstimator,
}

    #Obtain external states inputs for component
    #TODO: If converter has dynamics, need to reference states:
    #external_ix = device.input_port_mapping[device.converter]
    #vcvd = device_states[external_ix[1]]
    #vcvq = device_states[external_ix[2]]
    external_ix = get_input_port_ix(device, PSY.LCLFilter)
    δ = device_states[external_ix[1]]

    #Obtain inner variables for component
    V_tR = get_inner_vars(device)[VR_inv_var]
    V_tI = get_inner_vars(device)[VI_inv_var]
    vcvd = get_inner_vars(device)[Vd_cnv_var]
    vcvq = get_inner_vars(device)[Vq_cnv_var]

    #Get parameters
    filter = PSY.get_filter(device)
    ωb = 2 * pi * f0
    lf = PSY.get_lf(filter)
    rf = PSY.get_rf(filter)
    cf = PSY.get_cf(filter)
    lg = PSY.get_lg(filter)
    rg = PSY.get_rg(filter)
    MVABase = PSY.get_inverter_Sbase(device)

    #RI to dq transformation
    V_dq = ri_dq(δ) * [V_tR; V_tI]
    V_g = sqrt(V_tR^2 + V_tI^2)

    #Obtain indices for component w/r to device
    local_ix = get_local_state_ix(device, PSY.LCLFilter)

    #Define internal states for filter
    internal_states = @view device_states[local_ix]
    icvd = internal_states[1]
    icvq = internal_states[2]
    vod = internal_states[3]
    voq = internal_states[4]
    iod = internal_states[5]
    ioq = internal_states[6]

    #Inputs (control signals) - N/A

    #Compute 6 states ODEs (D'Arco EPSR122 Model)
    #Inverter Output Inductor (internal state)
    #𝜕id_c/𝜕t
    output_ode[local_ix[1]] =
        (ωb / lf * vcvd - ωb / lf * vod - ωb * rf / lf * icvd + ωb * ω_sys * icvq)
    #𝜕iq_c/𝜕t
    output_ode[local_ix[2]] =
        (ωb / lf * vcvq - ωb / lf * voq - ωb * rf / lf * icvq - ωb * ω_sys * icvd)
    #LCL Capacitor (internal state)
    #𝜕vd_o/𝜕t
    output_ode[local_ix[3]] = (ωb / cf * icvd - ωb / cf * iod + ωb * ω_sys * voq)
    #𝜕vq_o/𝜕t
    output_ode[local_ix[4]] = (ωb / cf * icvq - ωb / cf * ioq - ωb * ω_sys * vod)
    #Grid Inductance (internal state)
    #𝜕id_o/𝜕t
    output_ode[local_ix[5]] =
        (ωb / lg * vod - ωb / lg * V_dq[2] - ωb * rg / lg * iod + ωb * ω_sys * ioq)
    #𝜕iq_o/𝜕t
    output_ode[local_ix[6]] =
        (ωb / lg * voq + ωb / lg * V_dq[1] - ωb * rg / lg * ioq - ωb * ω_sys * iod)

    #Update inner_vars
    get_inner_vars(device)[Vd_filter_var] = vod
    get_inner_vars(device)[Vq_filter_var] = voq
    #TODO: If PLL models at PCC, need to update inner vars:
    #get_inner_vars(device)[Vd_filter_var] = V_dq[q::dq_ref]
    #get_inner_vars(device)[Vq_filter_var] = V_dq[q::dq_ref]

    #Compute current from the generator to the grid
    I_RI = (MVABase / sys_Sbase) * dq_ri(δ) * [iod; ioq]
    #Update current
    current_r[1] += I_RI[1]
    current_i[1] += I_RI[2]
end
