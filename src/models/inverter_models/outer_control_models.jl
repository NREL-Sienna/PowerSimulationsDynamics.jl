function mass_matrix_outer_entries!(
    mass_matrix,
    outer_control::O,
    global_index::Base.ImmutableDict{Symbol, Int64},
) where {O <: PSY.OuterControl}
    @debug "Using default mass matrix entries $O"
end

#####################################################
### Auxiliary ODE calculations via Flags dispatch ###
#####################################################

### Active Controllers ###

#Freq_Flag = 1
function _mdl_ode_RE_active_controller_AB!(
    active_controller_ode,
    active_controller_states,
    p_elec_out,
    ω_sys,
    Vt_filt,
    ::Type{Base.RefValue{1}},
    active_power_control::PSY.ActiveRenewableControllerAB,
    dynamic_device::DynamicWrapper{
        PSY.DynamicInverter{
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
    },
    inner_vars::AbstractVector,
) where {
    C <: PSY.Converter,
    IC <: PSY.InnerControl,
    DC <: PSY.DCSource,
    P <: PSY.FrequencyEstimator,
    F <: PSY.Filter,
}

    #Obtain external parameters
    p_ref = get_P_ref(dynamic_device)
    ω_ref = get_ω_ref(dynamic_device)
    # To do: Obtain proper frequency for a plant. For now using the system frequency.
    ω_plant = ω_sys

    #Obtain additional Active Power Controller parameters
    K_pg = PSY.get_K_pg(active_power_control)
    K_ig = PSY.get_K_ig(active_power_control)
    T_p = PSY.get_T_p(active_power_control)
    fdbd1, fdbd2 = PSY.get_fdbd_pnts(active_power_control)
    fe_min, fe_max = PSY.get_fe_lim(active_power_control)
    P_min, P_max = PSY.get_P_lim(active_power_control)
    T_g = PSY.get_T_g(active_power_control)
    D_dn = PSY.get_D_dn(active_power_control)
    D_up = PSY.get_D_up(active_power_control)
    dP_min, dP_max = PSY.get_dP_lim(active_power_control)
    P_min_inner, P_max_inner = PSY.get_P_lim_inner(active_power_control)
    T_pord = PSY.get_T_pord(active_power_control)

    #Define internal states for outer control
    p_flt = active_controller_states[1]
    ξ_P = active_controller_states[2]
    p_ext = active_controller_states[3]
    p_ord = active_controller_states[4]

    #Compute additional terms
    f_err = deadband_function(ω_ref - ω_plant, fdbd1, fdbd2)
    p_droop = min(D_dn * f_err, 0.0) + max(D_up * f_err, 0.0)
    p_err = clamp(p_ref + p_droop - p_flt, fe_min, fe_max)
    P_pi = K_pg * p_err + K_ig * ξ_P
    P_pi_binary = P_min <= P_pi <= P_max ? 1.0 : 0.0

    #Ramp and limits for p_ord
    #Hard limits
    p_ord_sat = clamp(p_ord, P_min_inner, P_max_inner)
    p_ord_binary = P_min_inner <= p_ord <= P_max_inner ? 1.0 : 0.0
    #Ramp limits: To check if 1/T_g has to be included or not
    p_ord_in = clamp(p_ext - p_ord_sat, dP_min, dP_max)

    #Update ODEs
    active_controller_ode[1] = (1.0 / T_p) * (p_elec_out - p_flt)
    active_controller_ode[2] = P_pi_binary * p_err
    active_controller_ode[3] = (1.0 / T_g) * (P_pi - p_ext)
    active_controller_ode[4] = p_ord_binary * (1.0 / T_pord) * p_ord_in

    #Update Inner Vars: Ioc_pcmd
    inner_vars[Id_oc_var] = p_ord_sat / max(Vt_filt, VOLTAGE_DIVISION_LOWER_BOUND)
end

#Freq_Flag = 0
function _mdl_ode_RE_active_controller_AB!(
    active_controller_ode,
    active_controller_states,
    p_elec_out,
    ω_sys,
    Vt_filt,
    ::Type{Base.RefValue{0}},
    active_power_control::PSY.ActiveRenewableControllerAB,
    dynamic_device::DynamicWrapper{
        PSY.DynamicInverter{
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
    },
    inner_vars::AbstractVector,
) where {
    C <: PSY.Converter,
    IC <: PSY.InnerControl,
    DC <: PSY.DCSource,
    P <: PSY.FrequencyEstimator,
    F <: PSY.Filter,
}

    #Obtain external parameters
    p_ref = get_P_ref(dynamic_device)
    #Obtain additional Active Power Controller parameters
    T_pord = PSY.get_T_pord(active_power_control)

    #Define internal states for outer control
    p_ord = active_controller_states[1]

    #Update ODE
    active_controller_ode[1] = (1.0 / T_pord) * (p_ref - p_ord)

    #Update Inner Vars: Ioc_pcmd
    inner_vars[Id_oc_var] = p_ord / max(Vt_filt, VOLTAGE_DIVISION_LOWER_BOUND)
end

### Reactive Controllers ###

# Order of Flags are:
# Ref_Flag, PF_Flag, V_Flag

#VC_Flag == N/A && Ref_Flag == 0 && PF_Flag == 0 && V_Flag == 1
#Named: Fixed Plant Q: Not compliant if Q_Flag = 0 from Inner Control
#Named: Plant Q and Local Q/V if Q_Flag = 1.
function _mdl_ode_RE_reactive_controller_AB!(
    reactive_controller_ode,
    reactive_controller_states,
    q_elec_out,
    Vt_filt,
    ::Type{Base.RefValue{0}},
    ::Type{Base.RefValue{0}},
    ::Type{Base.RefValue{1}},
    reactive_power_control::PSY.ReactiveRenewableControllerAB,
    dynamic_device::DynamicWrapper{
        PSY.DynamicInverter{
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
    },
    inner_vars::AbstractVector,
) where {
    C <: PSY.Converter,
    IC <: PSY.InnerControl,
    DC <: PSY.DCSource,
    P <: PSY.FrequencyEstimator,
    F <: PSY.Filter,
}

    #Obtain external parameters
    q_ref = get_Q_ref(dynamic_device)

    # Get Reactive Controller parameters
    T_fltr = PSY.get_T_fltr(reactive_power_control)
    K_p = PSY.get_K_p(reactive_power_control)
    K_i = PSY.get_K_i(reactive_power_control)
    T_ft = PSY.get_T_ft(reactive_power_control)
    T_fv = PSY.get_T_fv(reactive_power_control)
    # V_frz not implemented yet
    # V_frz = PSY.get_V_frz(reactive_power_control)
    e_min, e_max = PSY.get_e_lim(reactive_power_control)
    dbd1, dbd2 = PSY.get_dbd_pnts(reactive_power_control)
    Q_min, Q_max = PSY.get_Q_lim(reactive_power_control)
    Q_min_inner, Q_max_inner = PSY.get_Q_lim_inner(reactive_power_control)
    V_min, V_max = PSY.get_V_lim(reactive_power_control)
    K_qp = PSY.get_K_qp(reactive_power_control)
    K_qi = PSY.get_K_qi(reactive_power_control)

    #Define internal states for Reactive Control
    q_flt = reactive_controller_states[1]
    ξq_oc = reactive_controller_states[2]
    q_LL = reactive_controller_states[3]
    ξ_Q = reactive_controller_states[4]

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
    reactive_controller_ode[1] = (1.0 / T_fltr) * (q_elec_out - q_flt)
    reactive_controller_ode[2] = Q_binary_logic * q_err
    reactive_controller_ode[3] = (1.0 / T_fv) * (Q_pi_sat * (1.0 - T_ft / T_fv) - q_LL)
    reactive_controller_ode[4] = V_binary_logic * V_pi_in

    #Update Inner Vars
    inner_vars[V_oc_var] = V_pi_sat - Vt_filt
    inner_vars[Iq_oc_var] = Q_ext / max(Vt_filt, VOLTAGE_DIVISION_LOWER_BOUND)
end

#VC_Flag == N/A && Ref_Flag == 0 && PF_Flag == 0 && V_Flag == 0
#Fixed Plant Q if Q_Flag = 0 from Inner Control (not Compliant from CAISO requirements)
function _mdl_ode_RE_reactive_controller_AB!(
    reactive_controller_ode,
    reactive_controller_states,
    q_elec_out,
    Vt_filt,
    ::Type{Base.RefValue{0}},
    ::Type{Base.RefValue{0}},
    ::Type{Base.RefValue{0}},
    reactive_power_control::PSY.ReactiveRenewableControllerAB,
    dynamic_device::DynamicWrapper{
        PSY.DynamicInverter{
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
    },
    inner_vars::AbstractVector,
) where {
    C <: PSY.Converter,
    IC <: PSY.InnerControl,
    DC <: PSY.DCSource,
    P <: PSY.FrequencyEstimator,
    F <: PSY.Filter,
}
    #Obtain external parameters
    q_ref = get_Q_ref(dynamic_device)

    # Get Reactive Controller parameters
    T_fltr = PSY.get_T_fltr(reactive_power_control)
    K_p = PSY.get_K_p(reactive_power_control)
    K_i = PSY.get_K_i(reactive_power_control)
    T_ft = PSY.get_T_ft(reactive_power_control)
    T_fv = PSY.get_T_fv(reactive_power_control)
    # V_frz not implemented yet
    # V_frz = PSY.get_V_frz(reactive_power_control)
    e_min, e_max = PSY.get_e_lim(reactive_power_control)
    dbd1, dbd2 = PSY.get_dbd_pnts(reactive_power_control)
    Q_min, Q_max = PSY.get_Q_lim(reactive_power_control)

    #Define internal states for Reactive Control
    q_flt = reactive_controller_states[1]
    ξq_oc = reactive_controller_states[2]
    q_LL = reactive_controller_states[3]

    #Compute additional variables
    q_err = clamp(deadband_function(q_ref - q_flt, dbd1, dbd2), e_min, e_max)
    #Q error - PI controller
    Q_pi = K_p * q_err + K_i * ξq_oc
    Q_pi_sat = clamp(Q_pi, Q_min, Q_max)
    Q_binary_logic = Q_min <= Q_pi <= Q_max ? 1.0 : 0.0
    #Lead-Lag block
    Q_ext = q_LL + (T_ft / T_fv) * Q_pi_sat

    #Update ODEs
    reactive_controller_ode[1] = (1.0 / T_fltr) * (q_elec_out - q_flt)
    reactive_controller_ode[2] = Q_binary_logic * q_err
    reactive_controller_ode[3] = (1.0 / T_fv) * (Q_pi_sat * (1.0 - T_ft / T_fv) - q_LL)

    #Update Inner Vars
    inner_vars[V_oc_var] = Q_ext - Vt_filt
    inner_vars[Iq_oc_var] = Q_ext / max(Vt_filt, VOLTAGE_DIVISION_LOWER_BOUND)
end

############################################
### ODE calculations via device dispatch ###
############################################

function mdl_outer_ode!(
    device_states,
    output_ode,
    inner_vars,
    ω_sys,
    dynamic_device::DynamicWrapper{
        PSY.DynamicInverter{
            C,
            PSY.OuterControl{PSY.VirtualInertia, PSY.ReactivePowerDroop},
            IC,
            DC,
            P,
            F,
        },
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
    ω_pll = inner_vars[ω_freq_estimator_var]

    #Get Active Power Controller parameters
    outer_control = PSY.get_outer_control(dynamic_device)
    active_power_control = PSY.get_active_power(outer_control)
    Ta = PSY.get_Ta(active_power_control) #VSM Inertia constant
    kd = PSY.get_kd(active_power_control) #VSM damping constant
    kω = PSY.get_kω(active_power_control) #Frequency droop gain
    f0 = get_system_base_frequency(dynamic_device)
    ωb = 2 * pi * f0 #Rated angular frequency

    #Get Reactive Power Controller parameters
    reactive_power_control = PSY.get_reactive_power(outer_control)
    kq = PSY.get_kq(reactive_power_control) #Reactive power droop gain
    ωf = PSY.get_ωf(reactive_power_control) #Reactive power filter cutoff frequency

    #Obtain external parameters
    p_ref = get_P_ref(dynamic_device)
    ω_ref = get_ω_ref(dynamic_device)
    V_ref = get_V_ref(dynamic_device)
    q_ref = get_Q_ref(dynamic_device)

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
    inner_vars[θ_oc_var] = θ_oc
    inner_vars[ω_oc_var] = ω_oc
    inner_vars[V_oc_var] = V_ref + kq * (q_ref - qm)
end

function mdl_outer_ode!(
    device_states,
    output_ode,
    inner_vars,
    ω_sys,
    dynamic_device::DynamicWrapper{
        PSY.DynamicInverter{
            C,
            PSY.OuterControl{PSY.ActivePowerDroop, PSY.ReactivePowerDroop},
            IC,
            DC,
            P,
            F,
        },
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
    V_tR = inner_vars[Vr_inv_var]
    V_tI = inner_vars[Vi_inv_var]

    #Get Active Power Controller parameters
    outer_control = PSY.get_outer_control(dynamic_device)
    active_power_control = PSY.get_active_power(outer_control)
    Rp = PSY.get_Rp(active_power_control) #Droop Gain
    ωz = PSY.get_ωz(active_power_control) #Frequency cutoff frequency
    f0 = get_system_base_frequency(dynamic_device)
    ωb = 2 * pi * f0 #Rated angular frequency

    #Get Reactive Power Controller parameters
    reactive_power_control = PSY.get_reactive_power(outer_control)
    kq = PSY.get_kq(reactive_power_control) #Reactive power droop gain
    ωf = PSY.get_ωf(reactive_power_control) #Reactive power filter cutoff frequency

    #Obtain external parameters
    p_ref = get_P_ref(dynamic_device)
    ω_ref = get_ω_ref(dynamic_device)
    V_ref = get_V_ref(dynamic_device)
    q_ref = get_Q_ref(dynamic_device)

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
    inner_vars[θ_oc_var] = θ_oc
    inner_vars[ω_oc_var] = ω_oc
    inner_vars[V_oc_var] = V_ref + kq * (q_ref - qm)
end

function mdl_outer_ode!(
    device_states,
    output_ode,
    inner_vars,
    ω_sys,
    dynamic_device::DynamicWrapper{
        PSY.DynamicInverter{
            C,
            PSY.OuterControl{PSY.ActivePowerPI, PSY.ReactivePowerPI},
            IC,
            DC,
            P,
            F,
        },
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
    θ_pll = inner_vars[θ_freq_estimator_var]
    ω_pll = inner_vars[ω_freq_estimator_var]

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
    p_ref = get_P_ref(dynamic_device)
    q_ref = get_Q_ref(dynamic_device)

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
    inner_vars[θ_oc_var] = θ_pll
    inner_vars[ω_oc_var] = ω_pll
    inner_vars[Iq_oc_var] = Kp_p * (p_ref - p_oc) + Ki_p * σp_oc
    inner_vars[Id_oc_var] = Kp_q * (q_ref - q_oc) + Ki_q * σq_oc
end

function mdl_outer_ode!(
    device_states,
    output_ode,
    inner_vars,
    ω_sys,
    dynamic_device::DynamicWrapper{
        PSY.DynamicInverter{
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
        PSY.OuterControl{
            PSY.ActiveRenewableControllerAB,
            PSY.ReactiveRenewableControllerAB,
        },
    )
    Vt_filt = device_states[external_ix[1]]

    #Monitoring power from other branch not supported.
    V_R = inner_vars[Vr_inv_var]
    V_I = inner_vars[Vi_inv_var]
    I_R = inner_vars[Ir_inv_var]
    I_I = inner_vars[Ii_inv_var]
    p_elec_out = I_R * V_R + I_I * V_I
    q_elec_out = -I_I * V_R + I_R * V_I

    #Get Active Power Controller parameters
    outer_control = PSY.get_outer_control(dynamic_device)
    active_power_control = PSY.get_active_power(outer_control)
    #Note: Monitoring power from other branch not supported.
    Freq_Flag = PSY.get_Freq_Flag(active_power_control) #Frequency Flag

    #Get Reactive Power Controller parameters
    reactive_power_control = PSY.get_reactive_power(outer_control)
    #Note: Monitoring power from other branch not supported.
    Ref_Flag = PSY.get_Ref_Flag(reactive_power_control)
    PF_Flag = PSY.get_PF_Flag(reactive_power_control)
    V_Flag = PSY.get_V_Flag(reactive_power_control)

    #Obtain indices for component w/r to device
    local_ix = get_local_state_ix(
        dynamic_device,
        PSY.OuterControl{
            PSY.ActiveRenewableControllerAB,
            PSY.ReactiveRenewableControllerAB,
        },
    )
    internal_states = @view device_states[local_ix]
    internal_ode = @view output_ode[local_ix]
    active_n_states = PSY.get_n_states(active_power_control)
    reactive_n_states = PSY.get_n_states(reactive_power_control)
    active_ix_range = 1:active_n_states
    reactive_ix_range = (active_n_states + 1):(active_n_states + reactive_n_states)
    active_states = @view internal_states[active_ix_range]
    reactive_states = @view internal_states[reactive_ix_range]
    active_ode = @view internal_ode[active_ix_range]
    reactive_ode = @view internal_ode[reactive_ix_range]

    #Dispatch active power controller
    _mdl_ode_RE_active_controller_AB!(
        active_ode,
        active_states,
        p_elec_out,
        ω_sys,
        Vt_filt,
        Base.RefValue{Freq_Flag},
        active_power_control,
        dynamic_device,
        inner_vars,
    )

    #Dispatch reactive power controller
    _mdl_ode_RE_reactive_controller_AB!(
        reactive_ode,
        reactive_states,
        q_elec_out,
        Vt_filt,
        Base.RefValue{Ref_Flag},
        Base.RefValue{PF_Flag},
        Base.RefValue{V_Flag},
        reactive_power_control,
        dynamic_device,
        inner_vars,
    )
end
