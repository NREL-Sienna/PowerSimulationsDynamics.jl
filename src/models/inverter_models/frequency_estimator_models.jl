function mdl_freq_estimator_ode!(device_states,
        output_ode,
        f0,
        device::DynInverter{C,O,VC,DC,PLL,F}) where {C <: Converter,
                                                   O <: OuterControl,
                                                   VC<: VSControl,
                                                   DC<: DCSource,
                                                   F <: Filter}

    #Obtain external states inputs for component
    external_ix = device.input_port_mapping[device.freq_estimator]
    vod = device_states[external_ix[1]]
    voq = device_states[external_ix[2]]
    δθ_vsm = device_states[external_ix[3]]

    #Obtain inner variables for component
       #vod = device.inner_vars[Vdo_var]
       #voq = device.inner_vars[Vqo_var]
    #δθ_vsm = device.inner_vars[δdqRI_var]

    #Get parameters
      ω_lp = device.freq_estimator.ω_lp
    kp_pll = device.freq_estimator.kp_pll
    ki_pll = device.freq_estimator.ki_pll
        ωb = 2.0*pi*f0
        ωg = 1.0 #TODO:  create getter later

    #Obtain indices for component w/r to device
    local_ix = device.local_state_ix[device.freq_estimator]

    #Define internal states for frequency estimator
    internal_states = @view device_states[local_ix]
    vpll_d = internal_states[1]
    vpll_q = internal_states[2]
     ϵ_pll = internal_states[3]
    δθ_pll = internal_states[4]

    #Inputs (control signals)

    #Compute 6 states ODEs (D'Arco EPSR122 Model)
    #Output Voltage LPF (internal state)
    #𝜕vpll_d/𝜕t, D'Arco ESPR122 eqn. 12
    output_ode[local_ix[1]] = ( ω_lp*vod*cos(δθ_pll-δθ_vsm)
                              + ω_lp*voq*sin(δθ_pll-δθ_vsm)
                              - ω_lp*vpll_d)
    #𝜕vpll_q/𝜕t, D'Arco ESPR122 eqn. 12
    output_ode[local_ix[2]] = (- ω_lp*vod*sin(δθ_pll-δθ_vsm)
                               + ω_lp*voq*cos(δθ_pll-δθ_vsm)
                               - ω_lp*vpll_q)
    #PI Integrator (internal state)
    #𝜕dϵ_pll/𝜕t, D'Arco ESPR122 eqn. 13
    output_ode[local_ix[3]] = atan(vpll_q/vpll_d)
    #PLL Frequency Deviation (internal state)
    #𝜕δθ_pll/𝜕t, D'Arco ESPR122 eqn. 15
    output_ode[local_ix[4]] = ( ωb*kp_pll*atan(vpll_q/vpll_d)
                              + ωb*ki_pll*ϵ_pll)

    #Update inner_vars
    #PLL frequency, D'Arco EPSR122 eqn. 16
    device.inner_vars[ω_freq_estimator_var] = ( kp_pll*atan(vpll_q/vpll_d)
                                              + ki_pll*ϵ_pll
                                              + 1.0 )
        #TODO: replace 1.0 w/ grid getter, or not needed if PCC?
end
