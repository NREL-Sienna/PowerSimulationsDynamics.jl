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
    get_inner_vars(dynamic_device)[P_ES_var] = p_elec_out
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
    get_inner_vars(dynamic_device)[θ_oc_var] = θ0_oc
    get_inner_vars(dynamic_device)[ω_oc_var] = PSY.get_ω_ref(dynamic_device)
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
    get_inner_vars(dynamic_device)[P_ES_var] = p_elec_out
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
    get_inner_vars(dynamic_device)[θ_oc_var] = θ0_oc
    get_inner_vars(dynamic_device)[ω_oc_var] = PSY.get_ω_ref(dynamic_device)
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
    get_inner_vars(dynamic_device)[P_ES_var] = p_elec_out
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
    get_inner_vars(dynamic_device)[θ_oc_var] = θ0_oc
    get_inner_vars(dynamic_device)[ω_oc_var] = PSY.get_ω_ref(dynamic_device)
    get_inner_vars(dynamic_device)[Id_oc_var] = I_dq_cnv[d]
    get_inner_vars(dynamic_device)[Iq_oc_var] = I_dq_cnv[q]
    #Update Q_ref. Initialization assumes q_ref = q_elec_out from PF solution
    PSY.get_ext(dynamic_device)[CONTROL_REFS][P_ref_index] = p_elec_out
    PSY.set_P_ref!(PSY.get_active_power(PSY.get_outer_control(dynamic_device)), p_elec_out)
    PSY.get_ext(dynamic_device)[CONTROL_REFS][Q_ref_index] = q_elec_out
    PSY.set_Q_ref!(
        PSY.get_reactive_power(PSY.get_outer_control(dynamic_device)),
        q_elec_out,
    )
end

function initialize_outer!(
    device_states,
    static::PSY.StaticInjection,
    dynamic_device::PSY.DynamicInverter{
        C,
        PSY.OuterControl{
            PSY.ActiveRenewableControllerAB,
            PSY.ReactiveRenewableControllerAB,
        },
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
    function get_value_I(v::Float64)
        return v
    end
    function get_value_I(v::Int)
        return v
    end
    function get_value_I(v::ForwardDiff.Dual)
        return v.value
    end

    V_R = get_inner_vars(dynamic_device)[Vr_inv_var]
    V_I = get_inner_vars(dynamic_device)[Vi_inv_var]
    I_R = get_value_I(get_inner_vars(dynamic_device)[Ir_inv_var])
    I_I = get_value_I(get_inner_vars(dynamic_device)[Ii_inv_var])
    V_t = sqrt(V_R^2 + V_I^2)
    p_elec_out = I_R * V_R + I_I * V_I
    q_elec_out = -I_I * V_R + I_R * V_I
    q_ref = PSY.get_ext(dynamic_device)[PSID.CONTROL_REFS][PSID.Q_ref_index]

    #Get Outer Controller parameters
    outer_control = PSY.get_outer_control(dynamic_device)
    active_power_control = PSY.get_active_power(outer_control)
    Freq_Flag = PSY.get_Freq_Flag(active_power_control) #Frequency Flag

    #Set state counter for variable number of states due to flags
    state_ct = 1

    #Update inner_vars
    get_inner_vars(dynamic_device)[P_ES_var] = p_elec_out

    #Obtain indices for component w/r to device
    local_ix = get_local_state_ix(
        dynamic_device,
        PSY.OuterControl{
            PSY.ActiveRenewableControllerAB,
            PSY.ReactiveRenewableControllerAB,
        },
    )
    internal_states = @view device_states[local_ix]

    if Freq_Flag == 1
        #Obtain Parameters
        K_ig = PSY.get_K_ig(active_power_control)
        #Update States
        internal_states[state_ct] = p_elec_out
        internal_states[state_ct + 1] = p_elec_out / K_ig
        internal_states[state_ct + 2] = p_elec_out #p_ext
        internal_states[state_ct + 3] = p_elec_out #p_ord
        state_ct += 4
        #Update Inner Vars
        get_inner_vars(dynamic_device)[Id_oc_var] = p_elec_out / max(V_t, 0.01)
    else
        #Update States
        internal_states[state_ct] = p_elec_out
        #Update Inner Vars
        get_inner_vars(dynamic_device)[Id_oc_var] = p_elec_out / max(V_t, 0.01)
        state_ct += 1
    end

    reactive_power_control = PSY.get_reactive_power(outer_control)
    #Note: Monitoring power from other branch not supported.
    VC_Flag = PSY.get_VC_Flag(reactive_power_control)
    Ref_Flag = PSY.get_Ref_Flag(reactive_power_control)
    PF_Flag = PSY.get_PF_Flag(reactive_power_control)
    V_Flag = PSY.get_V_Flag(reactive_power_control)
    #Update references
    if VC_Flag == 0 && Ref_Flag == 0 && PF_Flag == 0 && V_Flag == 1
        #Get Reactive Controller Parameters
        K_i = PSY.get_K_i(reactive_power_control)
        K_qi = PSY.get_K_qi(reactive_power_control)
        #Update states
        internal_states[state_ct] = q_elec_out
        internal_states[state_ct + 1] = q_elec_out / K_i
        internal_states[state_ct + 2] = q_elec_out
        internal_states[state_ct + 3] = V_t / K_qi
        state_ct += 4
        #Update Inner Vars
        get_inner_vars(dynamic_device)[V_oc_var] = 0.0
        get_inner_vars(dynamic_device)[Iq_oc_var] = q_elec_out / max(V_t, 0.01)
    elseif VC_Flag == 0 && Ref_Flag == 0 && PF_Flag == 0 && V_Flag == 0
        K_i = PSY.get_K_i(reactive_power_control)
        #Update states
        internal_states[state_ct] = q_ref
        internal_states[state_ct + 1] = q_ref / K_i
        internal_states[state_ct + 2] = q_ref
        state_ct += 3
        #Update Inner Vars
        get_inner_vars(dynamic_device)[V_oc_var] = q_ref - V_t
        get_inner_vars(dynamic_device)[Iq_oc_var] = q_ref / max(V_t, 0.01)
    else
        error("Flags for Generic Renewable Model not supported yet")
    end
end
