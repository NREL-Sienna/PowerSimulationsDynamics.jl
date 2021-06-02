function mass_matrix_outer_entries!(
    mass_matrix,
   dynamic_device::PSY.DynamicInverter{C, O, IC, DC, P, F},
) where {
    C <: PSY.Converter,
    O <: PSY.OuterControl,
    IC <: PSY.InnerControl,
    DC <: PSY.DCSource,
    P <: PSY.FrequencyEstimator,
    F <: PSY.Filter,
}
    @debug "Using default mass matrix entries $O $(PSY.get_name(dynamic_device))"
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
    V_tR = get_inner_vars(dynamic_device)[VR_inv_var]
    V_tI = get_inner_vars(dynamic_device)[VI_inv_var]

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
