function mdl_VScontrol_ode!(
    device_states,
    output_ode,
    device::PSY.DynamicInverter{C, O, PSY.CombinedVIwithVZ, DC, P, F},
) where {
    C <: PSY.Converter,
    O <: PSY.OuterControl,
    DC <: PSY.DCSource,
    P <: PSY.FrequencyEstimator,
    F <: PSY.Filter,
}

    #Obtain external states inputs for component
    external_ix = get_input_port_ix(device, PSY.CombinedVIwithVZ)
    iod = device_states[external_ix[1]] #TODO: Should be var referemce?
    ioq = device_states[external_ix[2]] #TODO: Should be state reference?
    icvd = device_states[external_ix[3]]
    icvq = device_states[external_ix[4]]
    vod = device_states[external_ix[5]]
    voq = device_states[external_ix[6]]

    #Obtain inner variables for component
    #vod = get_inner_vars(device)[Vdo_var] #TODO: Should be state reference?
    #voq = get_inner_vars(device)[Vqo_var] #TODO: Should be state reference?
    ω_vsm = get_inner_vars(device)[ω_control_var]
    v_refr = get_inner_vars(device)[v_control_var]
    vdc = get_inner_vars(device)[Vdc_var]

    #Get Voltage Controller parameters
    vscontrol = PSY.get_vscontrol(device)
    filter = PSY.get_filter(device)
    kpv = PSY.get_kpv(vscontrol)
    kiv = PSY.get_kiv(vscontrol)
    kffi = PSY.get_kffi(vscontrol)
    cf = PSY.get_cf(filter)
    rv = PSY.get_rv(vscontrol)
    lv = PSY.get_lv(vscontrol)

    #Get Current Controller parameters
    kpc = PSY.get_kpc(vscontrol)
    kic = PSY.get_kic(vscontrol)
    kffv = PSY.get_kffv(vscontrol)
    lf = PSY.get_lf(filter)
    ωad = PSY.get_ωad(vscontrol)
    kad = PSY.get_kad(vscontrol)

    #Obtain indices for component w/r to device
    local_ix = get_local_state_ix(device, PSY.CombinedVIwithVZ)

    #Define internal states for frequency estimator
    internal_states = @view device_states[local_ix]
    ξ_d = internal_states[1]
    ξ_q = internal_states[2]
    γ_d = internal_states[3]
    γ_q = internal_states[4]
    ϕ_d = internal_states[5]
    ϕ_q = internal_states[6]

    #Inputs (control signals)

    ### Compute 6 states ODEs (D'Arco EPSR122 Model) ###
    ## SRF Voltage Control w/ Virtual Impedance ##
    #Virtual Impedance - Computation but not DAE
    #v_refr = V_ref + kq*(q_ref - qm)
    vod_ref = (v_refr - rv * iod + ω_vsm * lv * ioq)
    voq_ref = (-rv * ioq - ω_vsm * lv * iod)
    #Output Control Signal - Links to SRF Current Controller
    icvd_ref = (kpv * (vod_ref - vod) + kiv * ξ_d - cf * ω_vsm * voq + kffi * iod)

    icvq_ref = (kpv * (voq_ref - voq) + kiv * ξ_q + cf * ω_vsm * vod + kffi * ioq)
    #Voltage Control ODEs
    #PI Integrator (internal state)
    output_ode[local_ix[1]] = (vod_ref - vod)
    output_ode[local_ix[2]] = (voq_ref - voq)

    ## SRF Current Control ##
    #Active Damping
    #vad_d = kad*(vod-ϕ_d)
    #vad_q = kad*(voq-ϕ_q)
    #References for Converter Output Voltage
    vcvd_ref = (
        kpc * (icvd_ref - icvd) + kic * γ_d - ω_vsm * lf * icvq + kffv * vod -
            kad * (vod - ϕ_d)
    )
    vcvq_ref = (
        kpc * (icvq_ref - icvq) + kic * γ_q + ω_vsm * lf * icvd + kffv * voq -
            kad * (voq - ϕ_q)
    )
    #Modulation Commands to Converter
    #md = vcvd_ref/vdc
    #mq = vcvq_ref/vdc
    #Current Control ODEs
    #PI Integrator (internal state)
    output_ode[local_ix[3]] = icvd_ref - icvd
    output_ode[local_ix[4]] = icvq_ref - icvq
    #Active Damping LPF (internal state)
    output_ode[local_ix[5]] = ωad * vod - ωad * ϕ_d
    output_ode[local_ix[6]] = ωad * voq - ωad * ϕ_q

    #Update inner_vars
    get_inner_vars(device)[md_var] = vcvd_ref / vdc
    get_inner_vars(device)[mq_var] = vcvq_ref / vdc

end
