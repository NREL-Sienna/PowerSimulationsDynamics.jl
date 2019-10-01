function mdl_filter_ode!(device_states,
        output_ode,
        current_r,
        current_i,
        sys_Sbase,
        f0,
        device::DynInverter{C,O,VC,DC,P,LCLFilter}) where {C <: Converter,
                                                   O <: OuterControl,
                                                   VC<: VSControl,
                                                   DC<: DCSource,
                                                   P <: FrequencyEstimator}

    #Obtain external states inputs for component
    #TODO: If converter has dynamics, need to reference states:
    #external_ix = device.input_port_mapping[device.converter]
    #vcvd = device_states[external_ix[1]]
    #vcvq = device_states[external_ix[2]]
    external_ix = device.input_port_mapping[device.filter]
    δ = device_states[external_ix[1]]

    #Obtain inner variables for component
    V_tR = device.inner_vars[VR_inv_var]
    V_tI = device.inner_vars[VI_inv_var]
    vcvd = device.inner_vars[Vdcnv_var]
    vcvq = device.inner_vars[Vqcnv_var]

    #Get parameters
    ωb = 2*pi*f0
    lf = device.filter.lf
    rf = device.filter.rf
    cf = device.filter.cf
    lg = device.filter.lg
    rg = device.filter.rg
    MVABase = get_inverter_Sbase(device)
    ωg = 1.0 #TODO: create getter later

    #RI to dq transformation
    V_dq = ri_dq(δ)*[V_tR; V_tI]
    V_g = sqrt(V_tR^2 + V_tI^2)

    #Obtain indices for component w/r to device
    local_ix = device.local_state_ix[device.filter]
    #@show local_ix

    #Define internal states for filter
    internal_states = @view device_states[local_ix]
    icvd = internal_states[1]
    icvq = internal_states[2]
     vod = internal_states[3]
     voq = internal_states[4]
     iod = internal_states[5]
     ioq = internal_states[6]

    #Inputs (control signals) - N/A

    #Compute 6 states ODEs (D'Arco EPSR122 Model)
    #Inverter Output Inductor (internal state)
    #𝜕id_c/𝜕t
    output_ode[local_ix[1]] = ( ωb/lf*vcvd
                              - ωb/lf*vod
                              - ωb*rf/lf*icvd
                              + ωb*ωg*icvq )
    #𝜕iq_c/𝜕t
    output_ode[local_ix[2]] = ( ωb/lf*vcvq
                              - ωb/lf*voq
                              - ωb*rf/lf*icvq
                              - ωb*ωg*icvd )
    #LCL Capacitor (internal state)
    #𝜕vd_o/𝜕t
    output_ode[local_ix[3]] = ( ωb/cf*icvd
                              - ωb/cf*iod  #i_gd was specified; use equivalent i_od
                              + ωb*ωg*voq )
    #𝜕vq_o/𝜕t
    output_ode[local_ix[4]] = ( ωb/cf*icvq
                              - ωb/cf*ioq  #i_gq was specified; use equivalent i_oq
                              - ωb*ωg*vod )
    #Grid Inductance (internal state)
    #𝜕id_o/𝜕t
    output_ode[local_ix[5]] = ( ωb/lg*vod
                              - ωb/lg*V_dq[2] #vgd
                              - ωb*rg/lg*iod
                              + ωb*ωg*ioq )
    #𝜕iq_o/𝜕t
    output_ode[local_ix[6]] = ( ωb/lg*voq
                              + ωb/lg*V_dq[1] #vgq
                              - ωb*rg/lg*ioq
                              - ωb*ωg*iod )

    #Update inner_vars
    device.inner_vars[Vdo_var] = vod
    device.inner_vars[Vqo_var] = voq
    #TODO: If PLL models at PCC, need to update inner vars:
    #device.inner_vars[Vdo_var] = V_dq[q::dq_ref]
    #device.inner_vars[Vqo_var] = V_dq[q::dq_ref]

    #Compute current from the generator to the grid
    I_RI = (MVABase/sys_Sbase)*dq_ri(δ)*[iod; ioq]
    #@show MVABase
    #@show sys_Sbase
    #Update current
    current_r[1] += I_RI[1]
    current_i[1] += I_RI[2]
end
