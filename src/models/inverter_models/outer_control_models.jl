function mdl_outer_ode!(
    device_states,
    output_ode,
    f0,
    ω_sys,
    device::PSY.DynamicInverter{
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
        device,
        PSY.OuterControl{PSY.VirtualInertia, PSY.ReactivePowerDroop},
    )
    vpll_d = device_states[external_ix[1]]
    vpll_q = device_states[external_ix[2]]
    ϵ_pll = device_states[external_ix[3]]
    Vd_filter = device_states[external_ix[4]] #TODO: Should be inner reference after initialization
    Vq_filter = device_states[external_ix[5]] #TODO: Should be inner reference after initialization
    Id_filter = device_states[external_ix[6]]
    Iq_filter = device_states[external_ix[7]]

    #Obtain inner variables for component
    ω_pll = get_inner_vars(device)[ω_freq_estimator_var]

    #Get Active Power Controller parameters
    outer_control = PSY.get_outercontrol(device)
    active_power_control = PSY.get_active_power(outer_control)
    Ta = PSY.get_Ta(active_power_control) #VSM Inertia constant
    kd = PSY.get_kd(active_power_control) #VSM damping constant
    kω = PSY.get_kω(active_power_control) #Frequency droop gain
    ωb = PSY.get_ωb(active_power_control) #Rated angular frequency

    #Get Reactive Power Controller parameters
    reactive_power_control = PSY.get_reactive_power(outer_control)
    kq = PSY.get_kq(reactive_power_control) #Reactive power droop gain
    ωf = PSY.get_ωf(reactive_power_control) #Reactive power filter cutoff frequency

    #Obtain external parameters
    pll_control = PSY.get_freq_estimator(device)
    kp_pll = PSY.get_kp_pll(pll_control)
    ki_pll = PSY.get_ki_pll(pll_control)
    p_ref = PSY.get_ext(device)[CONTROL_REFS][P_ref_index]
    ω_ref = PSY.get_ext(device)[CONTROL_REFS][ω_ref_index]
    V_ref = PSY.get_ext(device)[CONTROL_REFS][V_ref_index]
    q_ref = PSY.get_ext(device)[CONTROL_REFS][Q_ref_index]

    #Obtain indices for component w/r to device
    local_ix = get_local_state_ix(
        device,
        PSY.OuterControl{PSY.VirtualInertia, PSY.ReactivePowerDroop},
    )

    #Define internal states for frequency estimator
    internal_states = @view device_states[local_ix]
    ω_oc = internal_states[1]
    θ_oc = internal_states[2]
    qm = internal_states[3]

    #Obtain additional expressions
    p_elec_out = Id_filter * Vd_filter + Iq_filter * Vq_filter

    #Compute 3 states ODEs
    output_ode[local_ix[1]] = (
        p_ref / Ta - p_elec_out / Ta - kd * (ω_oc - ω_pll) / Ta - kω * (ω_oc - ω_ref) / Ta
    )
    output_ode[local_ix[2]] = ωb * (ω_oc - ω_sys)
    output_ode[local_ix[3]] = (-ωf * Iq_filter * Vd_filter + ωf * Id_filter * Vq_filter - ωf * qm)

    #Update inner vars
    get_inner_vars(device)[θ_oc_var] = θ_oc
    get_inner_vars(device)[ω_oc_var] = ω_oc
    get_inner_vars(device)[V_oc_var] = V_ref + kq * (q_ref - qm)
end
 #output_ode[local_ix[1]] = (
    #    -Id_filter * Vd_filter / Ta - Iq_filter * Vq_filter / Ta +
    #    kd * kp_pll * atan(vpll_q / vpll_d) / Ta +
    #    kd * ki_pll * ϵ_pll / Ta - (kd + kω) * (ω_oc - ω_sys) / Ta +
    #    p_ref / Ta +
    #    kω * ω_ref / Ta - kω * ω_sys / Ta
    #)
