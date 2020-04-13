function mdl_freq_estimator_ode!(
    device_states,
    output_ode,
    f0,
    ω_sys,
    device::PSY.DynamicInverter{C, O, IC, DC, PSY.KauraPLL, F},
) where {
    C <: PSY.Converter,
    O <: PSY.OuterControl,
    IC <: PSY.InnerControl,
    DC <: PSY.DCSource,
    F <: PSY.Filter,
}

    #Obtain external states inputs for component
    external_ix = get_input_port_ix(device, PSY.KauraPLL)
    Vd_filter = device_states[external_ix[1]]
    Vq_filter = device_states[external_ix[2]]
    θ_oc = device_states[external_ix[3]]

    #Obtain inner variables for component
    #Vd_filter = device.inner_vars[Vd_filter_var]
    #Vq_filter = device.inner_vars[Vq_filter_var]
    #θ_oc = device.inner_vars[θ_oc_var]

    #Get parameters
    pll_control = PSY.get_freq_estimator(device)
    ω_lp = PSY.get_ω_lp(pll_control)
    kp_pll = PSY.get_kp_pll(pll_control)
    ki_pll = PSY.get_ki_pll(pll_control)
    ωb = 2.0 * pi * f0

    #Obtain indices for component w/r to device
    local_ix = get_local_state_ix(device, PSY.KauraPLL)

    #Define internal states for frequency estimator
    internal_states = @view device_states[local_ix]
    vpll_d = internal_states[1]
    vpll_q = internal_states[2]
    ϵ_pll = internal_states[3]
    θ_pll = internal_states[4]

    #Inputs (control signals)

    #Compute 6 states ODEs (D'Arco EPSR122 Model)
    #Output Voltage LPF (internal state)
    #𝜕vpll_d/𝜕t, D'Arco ESPR122 eqn. 12
    output_ode[local_ix[1]] = (
        ω_lp * Vd_filter * cos(θ_pll - θ_oc) + ω_lp * Vq_filter * sin(θ_pll - θ_oc) -
        ω_lp * vpll_d
    )
    #𝜕vpll_q/𝜕t, D'Arco ESPR122 eqn. 12
    output_ode[local_ix[2]] = (
        -ω_lp * Vd_filter * sin(θ_pll - θ_oc) + ω_lp * Vq_filter * cos(θ_pll - θ_oc) -
        ω_lp * vpll_q
    )
    #PI Integrator (internal state)
    #𝜕dϵ_pll/𝜕t, D'Arco ESPR122 eqn. 13
    output_ode[local_ix[3]] = atan(vpll_q / vpll_d)
    #PLL Frequency Deviation (internal state)
    #𝜕θ_pll/𝜕t, D'Arco ESPR122 eqn. 15
    output_ode[local_ix[4]] = (ωb * kp_pll * atan(vpll_q / vpll_d) + ωb * ki_pll * ϵ_pll)

    #Update inner_vars
    #PLL frequency, D'Arco EPSR122 eqn. 16
    get_inner_vars(device)[ω_freq_estimator_var] =
        (kp_pll * atan(vpll_q / vpll_d) + ki_pll * ϵ_pll + ω_sys)
end
