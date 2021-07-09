function mass_matrix_outer_entries!(
    mass_matrix,
    outer_control::O,
    global_index::Dict{Symbol, Int64},
) where {O <: PSY.OuterControl}
    @debug "Using default mass matrix entries $O"
end

function mdl_outer_ode!(
    device_states,
    output_ode,
    f0,
    ω_sys,
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

    #Obtain inner variables for component
    ω_pll = get_inner_vars(dynamic_device)[ω_freq_estimator_var]

    #Get Active Power Controller parameters
    outer_control = PSY.get_outer_control(dynamic_device)
    active_power_control = PSY.get_active_power(outer_control)
    Ta = PSY.get_Ta(active_power_control) #VSM Inertia constant
    kd = PSY.get_kd(active_power_control) #VSM damping constant
    kω = PSY.get_kω(active_power_control) #Frequency droop gain
    ωb = 2 * pi * f0 #Rated angular frequency

    #Get Reactive Power Controller parameters
    reactive_power_control = PSY.get_reactive_power(outer_control)
    kq = PSY.get_kq(reactive_power_control) #Reactive power droop gain
    ωf = PSY.get_ωf(reactive_power_control) #Reactive power filter cutoff frequency

    #Obtain external parameters
    p_ref = PSY.get_ext(dynamic_device)[CONTROL_REFS][P_ref_index]
    ω_ref = PSY.get_ext(dynamic_device)[CONTROL_REFS][ω_ref_index]
    V_ref = PSY.get_ext(dynamic_device)[CONTROL_REFS][V_ref_index]
    q_ref = PSY.get_ext(dynamic_device)[CONTROL_REFS][Q_ref_index]

    #Obtain indices for component w/r to device
    local_ix = get_local_state_ix(
        dynamic_device,
        PSY.OuterControl{PSY.VirtualInertia, PSY.ReactivePowerDroop},
    )

    #Define internal states for Outer Control
    internal_states = @view device_states[local_ix]
    ω_oc = internal_states[1]
    θ_oc = internal_states[2]
    qm = internal_states[3]

    #Obtain additional expressions
    p_elec_out = Ir_filter * Vr_filter + Ii_filter * Vi_filter
    q_elec_out = -Ii_filter * Vr_filter + Ir_filter * Vi_filter

    #Compute 3 states ODEs
    output_ode[local_ix[1]] =
        (p_ref / Ta - p_elec_out / Ta - kd * (ω_oc - ω_pll) / Ta - kω * (ω_oc - ω_ref) / Ta)
    output_ode[local_ix[2]] = ωb * (ω_oc - ω_sys)
    output_ode[local_ix[3]] = (ωf * (q_elec_out - qm))

    #Update inner vars
    get_inner_vars(dynamic_device)[θ_oc_var] = θ_oc
    get_inner_vars(dynamic_device)[ω_oc_var] = ω_oc
    get_inner_vars(dynamic_device)[V_oc_var] = V_ref + kq * (q_ref - qm)
end

function mdl_outer_ode!(
    device_states,
    output_ode,
    f0,
    ω_sys,
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

    #Obtain inner variables for component
    V_tR = get_inner_vars(dynamic_device)[VR_inv_var]
    V_tI = get_inner_vars(dynamic_device)[VI_inv_var]

    #Get Active Power Controller parameters
    outer_control = PSY.get_outer_control(dynamic_device)
    active_power_control = PSY.get_active_power(outer_control)
    Rp = PSY.get_Rp(active_power_control) #Droop Gain
    ωz = PSY.get_ωz(active_power_control) #Frequency cutoff frequency
    ωb = 2 * pi * f0 #Rated angular frequency

    #Get Reactive Power Controller parameters
    reactive_power_control = PSY.get_reactive_power(outer_control)
    kq = PSY.get_kq(reactive_power_control) #Reactive power droop gain
    ωf = PSY.get_ωf(reactive_power_control) #Reactive power filter cutoff frequency

    #Obtain external parameters
    p_ref = PSY.get_ext(dynamic_device)[PSID.CONTROL_REFS][PSID.P_ref_index]
    ω_ref = PSY.get_ext(dynamic_device)[PSID.CONTROL_REFS][PSID.ω_ref_index]
    V_ref = PSY.get_ext(dynamic_device)[PSID.CONTROL_REFS][PSID.V_ref_index]
    q_ref = PSY.get_ext(dynamic_device)[PSID.CONTROL_REFS][PSID.Q_ref_index]

    #Obtain indices for component w/r to device
    local_ix = get_local_state_ix(
        dynamic_device,
        PSY.OuterControl{PSY.ActivePowerDroop, PSY.ReactivePowerDroop},
    )

    #Define internal states for frequency estimator
    internal_states = @view device_states[local_ix]
    θ_oc = internal_states[1]
    pm = internal_states[2]
    qm = internal_states[3]

    #Obtain additional expressions
    p_elec_out = Ir_filter * Vr_filter + Ii_filter * Vi_filter
    q_elec_out = -Ii_filter * Vr_filter + Ir_filter * Vi_filter

    #Compute Frequency from Droop
    ω_oc = ω_ref + Rp * (p_ref - pm)

    #Compute 3 states ODEs
    output_ode[local_ix[1]] = ωb * (ω_oc - ω_sys)
    output_ode[local_ix[2]] = (ωz * (p_elec_out - pm))
    output_ode[local_ix[3]] = (ωf * (q_elec_out - qm))

    #Update inner vars
    get_inner_vars(dynamic_device)[θ_oc_var] = θ_oc
    get_inner_vars(dynamic_device)[ω_oc_var] = ω_oc
    get_inner_vars(dynamic_device)[V_oc_var] = V_ref + kq * (q_ref - qm)
end

function mdl_outer_ode!(
    device_states,
    output_ode,
    f0,
    ω_sys,
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

    #Obtain inner variables for component
    θ_pll = get_inner_vars(dynamic_device)[θ_freq_estimator_var]
    ω_pll = get_inner_vars(dynamic_device)[ω_freq_estimator_var]

    #Get Active Power Controller parameters
    outer_control = PSY.get_outer_control(dynamic_device)
    active_power_control = PSY.get_active_power(outer_control)
    Kp_p = PSY.get_Kp_p(active_power_control) #Proportional Gain
    Ki_p = PSY.get_Ki_p(active_power_control) #Integral Gain
    ωz = PSY.get_ωz(active_power_control) #Frequency cutoff frequency

    #Get Reactive Power Controller parameters
    reactive_power_control = PSY.get_reactive_power(outer_control)
    Kp_q = PSY.get_Kp_q(reactive_power_control) #Proportional Gain
    Ki_q = PSY.get_Ki_q(reactive_power_control) #Integral Gain
    ωf = PSY.get_ωf(reactive_power_control) #Reactive power filter cutoff frequency

    #Obtain external parameters
    p_ref = PSY.get_ext(dynamic_device)[PSID.CONTROL_REFS][PSID.P_ref_index]
    q_ref = PSY.get_ext(dynamic_device)[PSID.CONTROL_REFS][PSID.Q_ref_index]

    #Obtain indices for component w/r to device
    local_ix = get_local_state_ix(
        dynamic_device,
        PSY.OuterControl{PSY.ActivePowerPI, PSY.ReactivePowerPI},
    )

    #Define internal states for outer control
    internal_states = @view device_states[local_ix]
    σp_oc = internal_states[1]
    p_oc = internal_states[2]
    σq_oc = internal_states[3]
    q_oc = internal_states[4]

    #Obtain additional expressions
    p_elec_out = Ir_filter * Vr_filter + Ii_filter * Vi_filter
    q_elec_out = -Ii_filter * Vr_filter + Ir_filter * Vi_filter

    #Compute 4 states ODEs
    output_ode[local_ix[1]] = p_ref - p_oc
    output_ode[local_ix[2]] = ωz * (p_elec_out - p_oc)
    output_ode[local_ix[3]] = q_ref - q_oc
    output_ode[local_ix[4]] = ωf * (q_elec_out - q_oc)

    #Update inner vars
    get_inner_vars(dynamic_device)[θ_oc_var] = θ_pll
    get_inner_vars(dynamic_device)[ω_oc_var] = ω_pll
    get_inner_vars(dynamic_device)[Iq_oc_var] = Kp_p * (p_ref - p_oc) + Ki_p * σp_oc
    get_inner_vars(dynamic_device)[Id_oc_var] = Kp_q * (q_ref - q_oc) + Ki_q * σq_oc
end



function mdl_outer_ode!(
    device_states,
    output_ode,
    f0,
    ω_sys,
    dynamic_device::PSY.DynamicInverter{
        C,
        PSY.OuterControl{PSY.ActiveRenewableTypeAB, PSY.ReactiveRenewableTypeAB},
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
        PSY.OuterControl{PSY.ActiveRenewableTypeAB, PSY.ReactiveRenewableTypeAB},
    )
    Vt_filt = device_states[external_ix[1]]

    #Monitoring power from other branch not supported.
    V_R = get_inner_vars(dynamic_device)[VR_inv_var]
    V_I = get_inner_vars(dynamic_device)[VI_inv_var]
    I_R = get_inner_vars(dynamic_device)[Id_cnv_var]
    I_I = get_inner_vars(dynamic_device)[Iq_cnv_var]
    p_elec_out = I_R * V_R + I_I * V_I
    q_elec_out = -I_I * V_R + I_R * V_I

    #Obtain external parameters
    p_ref = PSY.get_ext(dynamic_device)[PSID.CONTROL_REFS][PSID.P_ref_index]
    ω_ref = PSY.get_ext(dynamic_device)[PSID.CONTROL_REFS][PSID.ω_ref_index]
    V_ref = PSY.get_ext(dynamic_device)[PSID.CONTROL_REFS][PSID.V_ref_index]
    q_ref = PSY.get_ext(dynamic_device)[PSID.CONTROL_REFS][PSID.Q_ref_index]
    # To do: Obtain proper frequency for a plant. For now using the system frequency.
    ω_plant = ω_sys

    #Set state counter for variable number of states due to flags
    state_ct = 1

    #Get Active Power Controller parameters
    outer_control = PSY.get_outer_control(dynamic_device)
    active_power_control = PSY.get_active_power(outer_control)
    #Note: Monitoring power from other branch not supported.
    Freq_Flag = PSY.get_Freq_Flag(active_power_control) #Frequency Flag
    p_ext = p_ref

    #Obtain indices for component w/r to device
    local_ix = get_local_state_ix(
        dynamic_device,
        PSY.OuterControl{PSY.ActiveRenewableTypeAB, PSY.ReactiveRenewableTypeAB},
        )
    internal_states = @view device_states[local_ix]

    if Freq_Flag == 1
        #Obtain additional Active Power Controller parameters
        K_pg = PSY.get_K_pg(active_power_control)
        K_pi = PSY.get_K_pi(active_power_control)
        T_p = PSY.get_T_p(active_power_control)
        fdbd1 = PSY.get_fdbd1(active_power_control)
        fdbd2 = PSY.get_fdbd2(active_power_control)
        fe_min, fe_max = PSY.get_fe_lim(active_power_control)
        P_min, P_max = PSY.get_P_lim(active_power_control)
        T_g = PSY.get_T_g(active_power_control)
        D_dn = PSY.get_D_dn(active_power_control)
        D_up = PSY.get_D_up(active_power_control)
        dP_min, dP_max = PSY.get_dP_lim(active_power_control)
        P_min_inner, P_max_inner = PSY.get_P_lim_inner(active_power_control)

        #Define internal states for outer control
        p_flt = internal_states[state_ct]
        ξ_P = internal_states[state_ct + 1]
        p_ext = internal_states[state_ct + 2]
        p_ord = internal_states[state_ct + 3]

        #Compute additional terms
        f_err = deadband_function(ω_ref - ω_plant, fdbd1, fdbd2)
        p_droop = min(D_dn * f_err, 0.0) + max(D_up * f_err, 0.0)
        p_err = clamp(p_ref + p_droop - p_flt, fe_min, fe_max)
        # To do: Limiters for PI block
        P_pi = K_pg * p_err + K_pi * ξ_P
        P_pi_binary = P_min <= P_pi <= P_max ? 1.0 : 0.0

        #Ramp and limits for p_ord
        #Hard limits
        p_ord_sat = clamp(p_ord, P_min_inner, P_max_inner)
        p_ord_binary = P_min_inner <= p_ord <= P_max_inner ? 1.0 : 0.0
        #Ramp limits: To check if 1/T_g has to be included or not
        p_ord_in = clamp(p_ext - p_ord_sat, dP_min, dP_max)

        #Update ODEs
        output_ode[local_ix[state_ct]] = (1.0 / T_p) * (p_elec_out - p_flt)
        output_ode[local_ix[state_ct + 1]] = P_pi_binary * p_err
        output_ode[local_ix[state_ct + 2]] = (1.0 / T_g) * (P_pi - p_ext)
        output_ode[local_ix[state_ct + 3]] = p_ord_binary * (1.0 / T_pord) * p_ord_in
        state_ct += 4

        #Update Inner Vars: Ioc_pcmd
        get_inner_vars(dynamic_device)[Id_oc_var] = p_ord_sat / max(Vt_filt, 0.01)
    else
        #p_ext is fixed so no limits are applied
        output_ode[local_ix[state_ct]] = (1.0 / T_pord) * (p_ext - p_ord)
        state_ct += 1
    end
    reactive_power_control = PSY.get_reactive_power(outer_control)
    #Note: Monitoring power from other branch not supported.
    VC_Flag = PSY.get_VC_Flag(reactive_power_control)
    Ref_Flag = PSY.get_Ref_Flag(reactive_power_control)
    PF_Flag = PSY.get_PF_Flag(reactive_power_control)
    V_Flag = PSY.get_V_Flag(reactive_power_control)
    if VC_Flag == 0 && Ref_Flag == 0 && PF_Flag == 0 && V_Flag == 1
        # Get Reactive Controller parameters
        T_fltr = PSY.get_T_fltr(reactive_power_control)
        K_p = PSY.get_K_p(reactive_power_control)
        K_i = PSY.get_K_i(reactive_power_control)
        T_ft = PSY.get_T_ft(reactive_power_control)
        T_fv = PSY.get_T_fv(reactive_power_control)
        # V_frz not implemented yet
        V_frz = PSY.get_V_frz(reactive_power_control)
        e_min, e_max = PSY.get_e_lim(reactive_power_control)
        dbd1 = PSY.get_dbd1(reactive_power_control)
        dbd2 = PSY.get_dbd2(reactive_power_control)
        Q_min, Q_max = PSY.get_Q_lim(reactive_power_control)
        T_p = PSY.get_T_p(reactive_power_control)
        Q_min_inner, Q_max_inner =  PSY.get_Q_lim_inner(reactive_power_control)
        V_min, V_max = PSY.get_V_lim(reactive_power_control)
        K_qp = PSY.get_K_qp(reactive_power_control)
        K_qi = PSY.get_K_qp(reactive_power_control)

        #Define internal states for Reactive Control
        q_flt = internal_states[state_ct]
        ξq_oc = internal_states[state_ct + 1]
        q_LL = internal_states[state_ct + 2]
        ξ_Q = internal_states[state_ct + 3]

        #Compute additional variables
        q_err = clamp(deadband_function(q_ref - q_flt, dbd1, dbd2), e_min, e_max)
        #Q error - PI controller
        Q_pi = K_p * q_err + K_i * ξq_oc
        Q_pi_sat = clamp(Q_pi, Q_min, Q_max)
        Q_binary_logic = Q_min <= Q_pi <= Q_max ? 1.0 : 0.0
        #Lead-Lag block
        Q_ext = q_LL + (T_ft / T_fv) * Q_pi_sat
        #PI voltage inner block
        V_pi_in = clamp(Q_ext, Q_min_inner, Q_max_inner) - q_elec_out
        V_pi = K_qp * V_pi_in + K_qi * ξ_Q
        V_pi_sat = clamp(V_pi, V_min, V_max)
        V_binary_logic = V_min <= V_pi <= V_max ? 1.0 : 0.0

        #Update ODEs
        output_ode[local_ix[state_ct]] = (1.0 / T_fltr) * (q_elec_out - q_flt)
        output_ode[local_ix[state_ct + 1]] = Q_binary_logic * q_err
        output_ode[local_ix[state_ct + 2]] = (1.0 / T_fv) * (Q_pi_sat * (1.0 - T_ft / T_fv) - q_LL) 
        output_ode[local_ix[state_ct + 3]] = V_binary_logic * V_pi_in
        state_ct += 4

        #Update Inner Vars
        get_inner_vars(dynamic_device)[V_oc_var] = V_pi_sat - Vt_filt
        get_inner_vars(dynamic_device)[Iq_oc_var] = Q_ext / max(Vt_filt, 0.01)
    else
        error("Flags for Generic Renewable Model not supported yet")
    end

    #Update inner vars
    get_inner_vars(dynamic_device)[θ_oc_var] = θ_pll
    get_inner_vars(dynamic_device)[ω_oc_var] = ω_pll
    get_inner_vars(dynamic_device)[Iq_oc_var] = Kp_p * (p_ref - p_oc) + Ki_p * σp_oc
    get_inner_vars(dynamic_device)[Id_oc_var] = Kp_q * (q_ref - q_oc) + Ki_q * σq_oc
end