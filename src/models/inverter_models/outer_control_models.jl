function mdl_outer_ode!(device_states,
                        output_ode,
                        f0,
                        device::PSY.DynamicInverter{C,PSY.VirtualInertiaQdroop{PSY.VirtualInertia,PSY.ReactivePowerDroop},VC,DC,P,F}) where {C <: PSY.Converter,
                                                                VC<: PSY.VSControl,
                                                                DC<: PSY.DCSource,
                                                                P <: PSY.FrequencyEstimator,
                                                                F <: PSY.Filter}


    #Obtain external states inputs for component
    external_ix = get_input_port_ix(device, PSY.VirtualInertiaQdroop{PSY.VirtualInertia,PSY.ReactivePowerDroop})
    vpll_d = device_states[external_ix[1]]
    vpll_q = device_states[external_ix[2]]
    ϵ_pll = device_states[external_ix[3]]
    vod = device_states[external_ix[4]]
    voq = device_states[external_ix[5]]
    iod = device_states[external_ix[6]]
    ioq = device_states[external_ix[7]]

    #Obtain inner variables for component
    ω_pll =  get_inner_vars(device)[ω_freq_estimator_var]

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
    p_ref = PSY.get_P_ref(device)
    ω_ref = PSY.get_ω_ref(device)
    V_ref = PSY.get_V_ref(device)
    #q_ref = PSY.get_Q_ref(device)
    q_ref = device.Q_ref
    ωg = 1.0

    #Obtain indices for component w/r to device
    local_ix = get_local_state_ix(device, PSY.VirtualInertiaQdroop{PSY.VirtualInertia,PSY.ReactivePowerDroop})

    #Define internal states for frequency estimator
    internal_states = @view device_states[local_ix]
    δω_vsm = internal_states[1]
    δθ_vsm = internal_states[2]
    qm = internal_states[3]

    #Compute 3 states ODEs
    output_ode[local_ix[1]] = (- iod*vod/Ta
                               - ioq*voq/Ta
                               + kd*kp_pll*atan(vpll_q/vpll_d)/Ta
                               + kd*ki_pll*ϵ_pll/Ta
                               - (kd+kω)*δω_vsm/Ta
                               + p_ref/Ta
                               + kω*ω_ref/Ta
                               - kω*ωg/Ta)
    output_ode[local_ix[2]] = ωb*δω_vsm
    output_ode[local_ix[3]] = (- ωf*ioq*vod
                               + ωf*iod*voq
                               - ωf*qm)

    #Update inner vars
    get_inner_vars(device)[δdqRI_var] = δθ_vsm
    get_inner_vars(device)[ω_control_var] = δω_vsm + 1.0
    get_inner_vars(device)[v_control_var] = V_ref + kq*(q_ref - qm)
end
