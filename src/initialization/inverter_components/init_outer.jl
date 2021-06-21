function initialize_outer!(
    device_states,
    static::PSY.StaticInjection,
    dynamic_device::PSY.DynamicInverter{
        C,
        PSY.OuterControl{PSY.VirtualInertia, PSY.ReactivePowerDroop},
        IC,
        DC,
        P,
        F,
    },
) where {
    C <: PSY.Converter,
    IC <: PSY.InnerControl,
    DC <: PSY.DCSource,
    P <: PSY.FrequencyEstimator,
    F <: PSY.Filter,
}

    #Obtain external states inputs for component
    external_ix = get_input_port_ix(
        dynamic_device,
        PSY.OuterControl{PSY.VirtualInertia, PSY.ReactivePowerDroop},
    )
    Vr_filter = device_states[external_ix[1]]
    Vi_filter = device_states[external_ix[2]]
    Ir_filter = device_states[external_ix[3]]
    Ii_filter = device_states[external_ix[4]]

    Vr_cnv = get_inner_vars(dynamic_device)[Vr_cnv_var]
    Vi_cnv = get_inner_vars(dynamic_device)[Vi_cnv_var]
    θ0_oc = angle(Vr_cnv + 1im * Vi_cnv)

    #Obtain additional expressions
    p_elec_out = Ir_filter * Vr_filter + Ii_filter * Vi_filter
    q_elec_out = -Ii_filter * Vr_filter + Ir_filter * Vi_filter

    #Update inner_vars
    set_inner_vars!(dynamic_device, P_ES_var, p_elec_out)
    #Update states
    outer_ix = get_local_state_ix(
        dynamic_device,
        PSY.OuterControl{PSY.VirtualInertia, PSY.ReactivePowerDroop},
    )
    outer_states = @view device_states[outer_ix]
    outer_states[1] = PSY.get_ω_ref(dynamic_device) #ω
    outer_states[2] = θ0_oc #θ_oc
    outer_states[3] = q_elec_out #qm

    #Update inner vars
    set_inner_vars!(dynamic_device, θ_oc_var, θ0_oc)
    set_inner_vars!(dynamic_device,ω_oc_var, PSY.get_ω_ref(dynamic_device))
    #Update Q_ref. Initialization assumes q_ref = q_elec_out of PF solution
    PSY.get_ext(dynamic_device)[CONTROL_REFS][P_ref_index] = p_elec_out
    PSY.set_P_ref!(PSY.get_active_power(PSY.get_outer_control(dynamic_device)), p_elec_out)
    PSY.get_ext(dynamic_device)[CONTROL_REFS][Q_ref_index] = q_elec_out
end

function initialize_outer!(
    device_states,
    static::PSY.StaticInjection,
    dynamic_device::PSY.DynamicInverter{
        C,
        PSY.OuterControl{PSY.ActivePowerDroop, PSY.ReactivePowerDroop},
        IC,
        DC,
        P,
        F,
    },
) where {
    C <: PSY.Converter,
    IC <: PSY.InnerControl,
    DC <: PSY.DCSource,
    P <: PSY.FrequencyEstimator,
    F <: PSY.Filter,
}

    #Obtain external states inputs for component
    external_ix = get_input_port_ix(
        dynamic_device,
        PSY.OuterControl{PSY.ActivePowerDroop, PSY.ReactivePowerDroop},
    )
    Vr_filter = device_states[external_ix[1]]
    Vi_filter = device_states[external_ix[2]]
    Ir_filter = device_states[external_ix[3]]
    Ii_filter = device_states[external_ix[4]]

    Vr_cnv = get_inner_vars(dynamic_device)[Vr_cnv_var]
    Vi_cnv = get_inner_vars(dynamic_device)[Vi_cnv_var]
    θ0_oc = angle(Vr_cnv + 1im * Vi_cnv)

    #Obtain additional expressions
    p_elec_out = Ir_filter * Vr_filter + Ii_filter * Vi_filter
    q_elec_out = -Ii_filter * Vr_filter + Ir_filter * Vi_filter

    #Update inner_vars
    set_inner_vars!(dynamic_device, P_ES_var, p_elec_out)
    #Update states
    outer_ix = get_local_state_ix(
        dynamic_device,
        PSY.OuterControl{PSY.ActivePowerDroop, PSY.ReactivePowerDroop},
    )
    outer_states = @view device_states[outer_ix]
    outer_states[1] = θ0_oc #θ_oc
    outer_states[2] = p_elec_out #pm
    outer_states[3] = q_elec_out #qm

    #Update inner vars
    set_inner_vars!(dynamic_device, θ_oc_var, θ0_oc)
    set_inner_vars!(dynamic_device, ω_oc_var, PSY.get_ω_ref(dynamic_device))
    #Update Q_ref. Initialization assumes q_ref = q_elec_out of PF solution
    PSY.get_ext(dynamic_device)[CONTROL_REFS][P_ref_index] = p_elec_out
    PSY.set_P_ref!(PSY.get_active_power(PSY.get_outer_control(dynamic_device)), p_elec_out)
    PSY.get_ext(dynamic_device)[CONTROL_REFS][Q_ref_index] = q_elec_out
end

function initialize_outer!(
    device_states,
    static::PSY.StaticInjection,
    dynamic_device::PSY.DynamicInverter{
        C,
        PSY.OuterControl{PSY.ActivePowerPI, PSY.ReactivePowerPI},
        IC,
        DC,
        P,
        F,
    },
) where {
    C <: PSY.Converter,
    IC <: PSY.InnerControl,
    DC <: PSY.DCSource,
    P <: PSY.FrequencyEstimator,
    F <: PSY.Filter,
}

    #Obtain external states inputs for component
    external_ix = get_input_port_ix(
        dynamic_device,
        PSY.OuterControl{PSY.ActivePowerPI, PSY.ReactivePowerPI},
    )
    Vr_filter = device_states[external_ix[1]]
    Vi_filter = device_states[external_ix[2]]
    Ir_filter = device_states[external_ix[3]]
    Ii_filter = device_states[external_ix[4]]
    Ir_cnv = device_states[external_ix[5]]
    Ii_cnv = device_states[external_ix[6]]

    #Obtain additional expressions
    θ0_oc = get_inner_vars(dynamic_device)[θ_freq_estimator_var]
    I_dq_cnv = ri_dq(θ0_oc + pi / 2) * [Ir_cnv; Ii_cnv]
    p_elec_out = Ir_filter * Vr_filter + Ii_filter * Vi_filter
    q_elec_out = -Ii_filter * Vr_filter + Ir_filter * Vi_filter

    #Get Outer Controller parameters
    outer_control = PSY.get_outer_control(dynamic_device)
    active_power_control = PSY.get_active_power(outer_control)
    Ki_p = PSY.get_Ki_p(active_power_control) #Integral Gain
    reactive_power_control = PSY.get_reactive_power(outer_control)
    Ki_q = PSY.get_Ki_q(reactive_power_control) #Integral Gain

    #Update inner_vars
    set_inner_vars!(dynamic_device, P_ES_var, p_elec_out)
    #Update states
    outer_ix = get_local_state_ix(
        dynamic_device,
        PSY.OuterControl{PSY.ActivePowerPI, PSY.ReactivePowerPI},
    )
    outer_states = @view device_states[outer_ix]
    outer_states[1] = I_dq_cnv[q] / Ki_p #σp_oc
    outer_states[2] = p_elec_out #p_oc
    outer_states[3] = I_dq_cnv[d] / Ki_q #σq_oc
    outer_states[4] = q_elec_out #q_oc

    #Update inner vars
    set_inner_vars!(dynamic_device, θ_oc_var, θ0_oc)
    set_inner_vars!(dynamic_device, ω_oc_var, PSY.get_ω_ref(dynamic_device))
    set_inner_vars!(dynamic_device, Id_oc_var, I_dq_cnv[d])
    set_inner_vars!(dynamic_device, Iq_oc_var, I_dq_cnv[q])
    #Update Q_ref. Initialization assumes q_ref = q_elec_out from PF solution
    PSY.get_ext(dynamic_device)[CONTROL_REFS][P_ref_index] = p_elec_out
    PSY.set_P_ref!(PSY.get_active_power(PSY.get_outer_control(dynamic_device)), p_elec_out)
    PSY.get_ext(dynamic_device)[CONTROL_REFS][Q_ref_index] = q_elec_out
    PSY.set_Q_ref!(
        PSY.get_reactive_power(PSY.get_outer_control(dynamic_device)),
        q_elec_out,
    )
end
