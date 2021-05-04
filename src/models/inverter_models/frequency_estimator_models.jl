function mdl_freq_estimator_ode!(
    device_states,
    output_ode,
    f0,
    ω_sys,
    dynamic_device::PSY.DynamicInverter{C, O, IC, DC, PSY.KauraPLL, F},
) where {
    C <: PSY.Converter,
    O <: PSY.OuterControl,
    IC <: PSY.InnerControl,
    DC <: PSY.DCSource,
    F <: PSY.Filter,
}

    #Obtain external states inputs for component
    external_ix = get_input_port_ix(dynamic_device, PSY.KauraPLL)
    Vr_filter = device_states[external_ix[1]]
    Vi_filter = device_states[external_ix[2]]

    #V_tR = get_inner_vars(dynamic_device)[VR_inv_var]
    #V_tI = get_inner_vars(dynamic_device)[VI_inv_var]

    #Get parameters
    pll_control = PSY.get_freq_estimator(dynamic_device)
    ω_lp = PSY.get_ω_lp(pll_control)
    kp_pll = PSY.get_kp_pll(pll_control)
    ki_pll = PSY.get_ki_pll(pll_control)
    ωb = 2.0 * pi * f0

    #Obtain indices for component w/r to device
    local_ix = get_local_state_ix(dynamic_device, PSY.KauraPLL)

    #Define internal states for frequency estimator
    internal_states = @view device_states[local_ix]
    vpll_d = internal_states[1]
    vpll_q = internal_states[2]
    ϵ_pll = internal_states[3]
    θ_pll = internal_states[4]

    #Transform to internal dq-PLL reference frame
    V_dq_pll = ri_dq(θ_pll + pi / 2) * [Vr_filter; Vi_filter]

    #Inputs (control signals)

    #Compute 6 states ODEs (D'Arco EPSR122 Model)
    #Output Voltage LPF (internal state)
    #𝜕vpll_d/𝜕t, D'Arco ESPR122 eqn. 12
    output_ode[local_ix[1]] = ω_lp * (V_dq_pll[d] - vpll_d)
    #𝜕vpll_q/𝜕t, D'Arco ESPR122 eqn. 12
    output_ode[local_ix[2]] = ω_lp * (V_dq_pll[q] - vpll_q)
    #PI Integrator (internal state)
    #𝜕dϵ_pll/𝜕t, D'Arco ESPR122 eqn. 13
    output_ode[local_ix[3]] = atan(vpll_q / vpll_d)
    #PLL Frequency Deviation (internal state)
    #𝜕θ_pll/𝜕t, D'Arco ESPR122 eqn. 15
    output_ode[local_ix[4]] = (ωb * kp_pll * atan(vpll_q / vpll_d) + ωb * ki_pll * ϵ_pll)

    #Update inner_vars
    #PLL frequency, D'Arco EPSR122 eqn. 16
    get_inner_vars(dynamic_device)[ω_freq_estimator_var] =
        (kp_pll * atan(vpll_q / vpll_d) + ki_pll * ϵ_pll + ω_sys)
end


function mdl_freq_estimator_ode!(
    device_states,
    output_ode,
    f0,
    ω_sys,
    dynamic_device::PSY.DynamicInverter{C, O, IC, DC, PSY.FixedFrequency, F},
) where {
    C <: PSY.Converter,
    O <: PSY.OuterControl,
    IC <: PSY.InnerControl,
    DC <: PSY.DCSource,
    F <: PSY.Filter,
}

    #Get parameters
    pll_control = PSY.get_freq_estimator(dynamic_device)
    frequency = PSY.get_frequency(pll_control)

    #Update inner_vars
    #PLL frequency
    get_inner_vars(dynamic_device)[ω_freq_estimator_var] = frequency
end