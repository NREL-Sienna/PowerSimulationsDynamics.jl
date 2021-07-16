function mass_matrix_inner_entries!(
    mass_matrix,
    inner_control::IC,
    global_index::Base.ImmutableDict{Symbol, Int64},
) where {IC <: PSY.InnerControl}
    @debug "Using default mass matrix entries $IC"
end

function mdl_inner_ode!(
    device_states,
    output_ode,
    dynamic_device::PSY.DynamicInverter{C, O, PSY.VoltageModeControl, DC, P, F},
) where {
    C <: PSY.Converter,
    O <: PSY.OuterControl,
    DC <: PSY.DCSource,
    P <: PSY.FrequencyEstimator,
    F <: PSY.Filter,
}

    #Obtain external states inputs for component
    external_ix = get_input_port_ix(dynamic_device, PSY.VoltageModeControl)
    Ir_filter = device_states[external_ix[1]]
    Ii_filter = device_states[external_ix[2]]
    Ir_cnv = device_states[external_ix[3]]
    Ii_cnv = device_states[external_ix[4]]
    Vr_filter = device_states[external_ix[5]] #TODO: Should be inner reference after initialization
    Vi_filter = device_states[external_ix[6]] #TODO: Should be inner reference after initialization

    #Obtain inner variables for component
    ω_oc = get_inner_vars(dynamic_device)[ω_oc_var]
    θ_oc = get_inner_vars(dynamic_device)[θ_oc_var]
    v_refr = get_inner_vars(dynamic_device)[V_oc_var]
    Vdc = get_inner_vars(dynamic_device)[Vdc_var]

    #Get Voltage Controller parameters
    inner_control = PSY.get_inner_control(dynamic_device)
    filter = PSY.get_filter(dynamic_device)
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
    local_ix = get_local_state_ix(dynamic_device, PSY.VoltageModeControl)

    #Define internal states for frequency estimator
    internal_states = @view device_states[local_ix]
    ξ_d = internal_states[1]
    ξ_q = internal_states[2]
    γ_d = internal_states[3]
    γ_q = internal_states[4]
    ϕ_d = internal_states[5]
    ϕ_q = internal_states[6]

    #Transformations to dq frame
    I_dq_filter = ri_dq(θ_oc + pi / 2) * [Ir_filter; Ii_filter]
    I_dq_cnv = ri_dq(θ_oc + pi / 2) * [Ir_cnv; Ii_cnv]
    V_dq_filter = ri_dq(θ_oc + pi / 2) * [Vr_filter; Vi_filter]

    ### Compute 6 states ODEs (D'Arco EPSR122 Model) ###
    ## SRF Voltage Control w/ Virtual Impedance ##
    #Virtual Impedance
    Vd_filter_ref = (v_refr - rv * I_dq_filter[d] + ω_oc * lv * I_dq_filter[q])
    Vq_filter_ref = (-rv * I_dq_filter[q] - ω_oc * lv * I_dq_filter[d])

    #Voltage Control ODEs
    #PI Integrator (internal state)
    output_ode[local_ix[1]] = (Vd_filter_ref - V_dq_filter[d])
    output_ode[local_ix[2]] = (Vq_filter_ref - V_dq_filter[q])

    #Output Control Signal - Links to SRF Current Controller
    Id_cnv_ref = (
        kpv * (Vd_filter_ref - V_dq_filter[d]) + kiv * ξ_d - cf * ω_oc * V_dq_filter[q] + kffi * I_dq_filter[d]
    )

    Iq_cnv_ref = (
        kpv * (Vq_filter_ref - V_dq_filter[q]) +
        kiv * ξ_q +
        cf * ω_oc * V_dq_filter[d] +
        kffi * I_dq_filter[q]
    )

    ## SRF Current Control ##
    #Current Control ODEs
    #PI Integrator (internal state)
    output_ode[local_ix[3]] = Id_cnv_ref - I_dq_cnv[d]
    output_ode[local_ix[4]] = Iq_cnv_ref - I_dq_cnv[q]

    #References for Converter Output Voltage
    Vd_cnv_ref = (
        kpc * (Id_cnv_ref - I_dq_cnv[d]) + kic * γ_d - ω_oc * lf * I_dq_cnv[q] +
        kffv * V_dq_filter[d] - kad * (V_dq_filter[d] - ϕ_d)
    )
    Vq_cnv_ref = (
        kpc * (Iq_cnv_ref - I_dq_cnv[q]) +
        kic * γ_q +
        ω_oc * lf * I_dq_cnv[d] +
        kffv * V_dq_filter[q] - kad * (V_dq_filter[q] - ϕ_q)
    )

    #Active Damping LPF (internal state)
    output_ode[local_ix[5]] = ωad * V_dq_filter[d] - ωad * ϕ_d
    output_ode[local_ix[6]] = ωad * V_dq_filter[q] - ωad * ϕ_q

    #Update inner_vars
    #Modulation Commands to Converter
    get_inner_vars(dynamic_device)[md_var] = Vd_cnv_ref / Vdc
    get_inner_vars(dynamic_device)[mq_var] = Vq_cnv_ref / Vdc
end

function mdl_inner_ode!(
    device_states,
    output_ode,
    dynamic_device::PSY.DynamicInverter{C, O, PSY.CurrentModeControl, DC, P, F},
) where {
    C <: PSY.Converter,
    O <: PSY.OuterControl,
    DC <: PSY.DCSource,
    P <: PSY.FrequencyEstimator,
    F <: PSY.Filter,
}

    #Obtain external states inputs for component
    external_ix = get_input_port_ix(dynamic_device, PSY.CurrentModeControl)
    # Ir_filter = device_states[external_ix[1]]
    # Ii_filter = device_states[external_ix[2]]
    Ir_cnv = device_states[external_ix[3]]
    Ii_cnv = device_states[external_ix[4]]
    Vr_filter = device_states[external_ix[5]] #TODO: Should be inner reference after initialization
    Vi_filter = device_states[external_ix[6]] #TODO: Should be inner reference after initialization

    #Obtain inner variables for component
    ω_oc = get_inner_vars(dynamic_device)[ω_oc_var]
    θ_oc = get_inner_vars(dynamic_device)[θ_oc_var]
    Vdc = get_inner_vars(dynamic_device)[Vdc_var]

    #Get Current Controller parameters
    inner_control = PSY.get_inner_control(dynamic_device)
    filter = PSY.get_filter(dynamic_device)
    kpc = PSY.get_kpc(inner_control)
    kic = PSY.get_kic(inner_control)
    kffv = PSY.get_kffv(inner_control)
    lf = PSY.get_lf(filter)

    #Obtain indices for component w/r to device
    local_ix = get_local_state_ix(dynamic_device, PSY.CurrentModeControl)

    #Define internal states for frequency estimator
    internal_states = @view device_states[local_ix]
    γ_d = internal_states[1]
    γ_q = internal_states[2]

    #Transformations to dq frame
    I_dq_cnv = ri_dq(θ_oc + pi / 2) * [Ir_cnv; Ii_cnv]
    V_dq_filter = ri_dq(θ_oc + pi / 2) * [Vr_filter; Vi_filter]

    #Input Control Signal - Links to SRF Current Controller
    Id_cnv_ref = get_inner_vars(dynamic_device)[Id_oc_var]
    Iq_cnv_ref = get_inner_vars(dynamic_device)[Iq_oc_var]

    #Current Control ODEs
    #PI Integrator (internal state)
    output_ode[local_ix[1]] = Id_cnv_ref - I_dq_cnv[d]
    output_ode[local_ix[2]] = Iq_cnv_ref - I_dq_cnv[q]

    #References for Converter Output Voltage
    Vd_cnv_ref = (
        kpc * (Id_cnv_ref - I_dq_cnv[d]) + kic * γ_d - ω_oc * lf * I_dq_cnv[q] +
        kffv * V_dq_filter[d]
    )
    Vq_cnv_ref = (
        kpc * (Iq_cnv_ref - I_dq_cnv[q]) +
        kic * γ_q +
        ω_oc * lf * I_dq_cnv[d] +
        kffv * V_dq_filter[q]
    )

    #Update inner_vars
    #Modulation Commands to Converter
    get_inner_vars(dynamic_device)[md_var] = Vd_cnv_ref / Vdc
    get_inner_vars(dynamic_device)[mq_var] = Vq_cnv_ref / Vdc
end
