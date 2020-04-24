function mdl_inner_ode!(
    device_states,
    output_ode,
    device::PSY.DynamicInverter{C, O, PSY.CurrentControl, DC, P, F},
) where {
    C <: PSY.Converter,
    O <: PSY.OuterControl,
    DC <: PSY.DCSource,
    P <: PSY.FrequencyEstimator,
    F <: PSY.Filter,
}

    #Obtain external states inputs for component
    external_ix = get_input_port_ix(device, PSY.CurrentControl)
    Ir_filter = device_states[external_ix[1]]
    Ii_filter = device_states[external_ix[2]]
    Ir_cnv = device_states[external_ix[3]]
    Ii_cnv = device_states[external_ix[4]]
    Vr_filter = device_states[external_ix[5]] #TODO: Should be inner reference after initialization
    Vi_filter = device_states[external_ix[6]] #TODO: Should be inner reference after initialization

    #Obtain inner variables for component
    #Vd_filter = get_inner_vars(device)[Vd_filter_var]
    #Vq_filter = get_inner_vars(device)[Vq_filter_var]
    ω_oc = get_inner_vars(device)[ω_oc_var]
    θ_oc = get_inner_vars(device)[θ_oc_var]
    v_refr = get_inner_vars(device)[V_oc_var]
    vdc = get_inner_vars(device)[Vdc_var]

    #Get Voltage Controller parameters
    inner_control = PSY.get_inner_control(device)
    filter = PSY.get_filter(device)
    kpv = PSY.get_kpv(inner_control)
    kiv = PSY.get_kiv(inner_control)
    kffi = PSY.get_kffi(inner_control)
    cf = PSY.get_cf(filter)
    rv = PSY.get_rv(inner_control)
    lv = PSY.get_lv(inner_control)

    #Get Current Controller parameters
    kpc = PSY.get_kpc(inner_control)
    kic = PSY.get_kic(inner_control)
    kffv = PSY.get_kffv(inner_control)
    lf = PSY.get_lf(filter)
    ωad = PSY.get_ωad(inner_control)
    kad = PSY.get_kad(inner_control)

    #Obtain indices for component w/r to device
    local_ix = get_local_state_ix(device, PSY.CurrentControl)

    #Define internal states for frequency estimator
    internal_states = @view device_states[local_ix]
    ξ_d = internal_states[1]
    ξ_q = internal_states[2]
    γ_d = internal_states[3]
    γ_q = internal_states[4]
    ϕ_d = internal_states[5]
    ϕ_q = internal_states[6]

    #Transformations
    I_dq_filter = ri_dq(θ_oc + pi / 2) * [Ir_filter; Ii_filter]
    I_dq_cnv = ri_dq(θ_oc + pi / 2) * [Ir_cnv; Ii_cnv]
    V_dq_filter = ri_dq(θ_oc + pi / 2) * [Vr_filter; Vi_filter]
    Id_filter = I_dq_filter[1]
    Iq_filter = I_dq_filter[2]
    Id_cnv = I_dq_cnv[1]
    Iq_cnv = I_dq_cnv[2]
    Vd_filter = V_dq_filter[1]
    Vq_filter = V_dq_filter[2]

    #Inputs (control signals)

    ### Compute 6 states ODEs (D'Arco EPSR122 Model) ###
    ## SRF Voltage Control w/ Virtual Impedance ##
    #Virtual Impedance - Computation but not DAE
    #v_refr = V_ref + kq*(q_ref - qm)
    Vd_filter_ref = (v_refr - rv * Id_filter + ω_oc * lv * Iq_filter)
    Vq_filter_ref = (-rv * Iq_filter - ω_oc * lv * Id_filter)
    #Output Control Signal - Links to SRF Current Controller
    Id_cnv_ref = (
        kpv * (Vd_filter_ref - Vd_filter) + kiv * ξ_d - cf * ω_oc * Vq_filter +
        kffi * Id_filter
    )

    Iq_cnv_ref = (
        kpv * (Vq_filter_ref - Vq_filter) +
        kiv * ξ_q +
        cf * ω_oc * Vd_filter +
        kffi * Iq_filter
    )
    #Voltage Control ODEs
    #PI Integrator (internal state)
    output_ode[local_ix[1]] = (Vd_filter_ref - Vd_filter)
    output_ode[local_ix[2]] = (Vq_filter_ref - Vq_filter)

    ## SRF Current Control ##
    #Active Damping
    #vad_d = kad*(Vd_filter-ϕ_d)
    #vad_q = kad*(Vq_filter-ϕ_q)
    #References for Converter Output Voltage
    Vd_cnv_ref = (
        kpc * (Id_cnv_ref - Id_cnv) + kic * γ_d - ω_oc * lf * Iq_cnv + kffv * Vd_filter - kad * (Vd_filter - ϕ_d)
    )
    Vq_cnv_ref = (
        kpc * (Iq_cnv_ref - Iq_cnv) + kic * γ_q + ω_oc * lf * Id_cnv + kffv * Vq_filter - kad * (Vq_filter - ϕ_q)
    )
    #Modulation Commands to Converter
    #md = Vd_cnv_ref/vdc
    #mq = Vq_cnv_ref/vdc
    #Current Control ODEs
    #PI Integrator (internal state)
    output_ode[local_ix[3]] = Id_cnv_ref - Id_cnv
    output_ode[local_ix[4]] = Iq_cnv_ref - Iq_cnv
    #Active Damping LPF (internal state)
    output_ode[local_ix[5]] = ωad * Vd_filter - ωad * ϕ_d
    output_ode[local_ix[6]] = ωad * Vq_filter - ωad * ϕ_q

    #Update inner_vars
    get_inner_vars(device)[md_var] = Vd_cnv_ref / vdc
    get_inner_vars(device)[mq_var] = Vq_cnv_ref / vdc

end
