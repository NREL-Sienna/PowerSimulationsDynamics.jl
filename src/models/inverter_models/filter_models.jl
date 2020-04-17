function mdl_filter_ode!(
    device_states,
    output_ode,
    current_r,
    current_i,
    sys_Sbase,
    f0,
    ω_sys,
    device::PSY.DynamicInverter{C, O, IC, DC, P, PSY.LCLFilter},
) where {
    C <: PSY.Converter,
    O <: PSY.OuterControl,
    IC <: PSY.InnerControl,
    DC <: PSY.DCSource,
    P <: PSY.FrequencyEstimator,
}

    #Obtain external states inputs for component
    #TODO: If converter has dynamics, need to reference states:
    #external_ix = device.input_port_mapping[device.converter]
    #Vd_cnv = device_states[external_ix[1]]
    #Vq_cnv = device_states[external_ix[2]]
    external_ix = get_input_port_ix(device, PSY.LCLFilter)
    δ = device_states[external_ix[1]]

    #Obtain inner variables for component
    V_tR = get_inner_vars(device)[VR_inv_var]
    V_tI = get_inner_vars(device)[VI_inv_var]
    Vd_cnv = get_inner_vars(device)[Vd_cnv_var]
    Vq_cnv = get_inner_vars(device)[Vq_cnv_var]

    #Get parameters
    filter = PSY.get_filter(device)
    ωb = 2 * pi * f0
    lf = PSY.get_lf(filter)
    rf = PSY.get_rf(filter)
    cf = PSY.get_cf(filter)
    lg = PSY.get_lg(filter)
    rg = PSY.get_rg(filter)
    MVABase = PSY.get_inverter_Sbase(device)

    #RI to dq transformation
    V_dq = ri_dq(δ+pi/2) * [V_tR; V_tI]
    V_g = sqrt(V_tR^2 + V_tI^2)

    #Obtain indices for component w/r to device
    local_ix = get_local_state_ix(device, PSY.LCLFilter)

    #Define internal states for filter
    internal_states = @view device_states[local_ix]
    Id_cnv = internal_states[1]
    Iq_cnv = internal_states[2]
    Vd_filter = internal_states[3]
    Vq_filter = internal_states[4]
    Id_filter = internal_states[5]
    Iq_filter = internal_states[6]

    #Inputs (control signals) - N/A

    #Compute 6 states ODEs (D'Arco EPSR122 Model)
    #Inverter Output Inductor (internal state)
    #𝜕id_c/𝜕t
    output_ode[local_ix[1]] = (
        ωb / lf * Vd_cnv - ωb / lf * Vd_filter - ωb * rf / lf * Id_cnv +
        ωb * ω_sys * Iq_cnv
    )
    #𝜕iq_c/𝜕t
    output_ode[local_ix[2]] = (
        ωb / lf * Vq_cnv - ωb / lf * Vq_filter - ωb * rf / lf * Iq_cnv -
        ωb * ω_sys * Id_cnv
    )
    #LCL Capacitor (internal state)
    #𝜕vd_o/𝜕t
    output_ode[local_ix[3]] =
        (ωb / cf * Id_cnv - ωb / cf * Id_filter + ωb * ω_sys * Vq_filter)
    #𝜕vq_o/𝜕t
    output_ode[local_ix[4]] =
        (ωb / cf * Iq_cnv - ωb / cf * Iq_filter - ωb * ω_sys * Vd_filter)
    #Grid Inductance (internal state)
    #𝜕id_o/𝜕t
    output_ode[local_ix[5]] = (
        ωb / lg * Vd_filter - ωb / lg * V_dq[1] - ωb * rg / lg * Id_filter +
        ωb * ω_sys * Iq_filter
    )
    #𝜕iq_o/𝜕t (Multiply Vq by -1 to lag instead of lead)
    output_ode[local_ix[6]] = (
        ωb / lg * Vq_filter + ωb / lg * (-V_dq[2]) - ωb * rg / lg * Iq_filter -
        ωb * ω_sys * Id_filter
    )

    #Update inner_vars
    get_inner_vars(device)[Vd_filter_var] = Vd_filter
    get_inner_vars(device)[Vq_filter_var] = Vq_filter
    #TODO: If PLL models at PCC, need to update inner vars:
    #get_inner_vars(device)[Vd_filter_var] = V_dq[q::dq_ref]
    #get_inner_vars(device)[Vq_filter_var] = V_dq[q::dq_ref]

    #Compute current from the inverter to the grid
    I_RI = (MVABase / sys_Sbase) * dq_ri(δ+pi/2) * [Id_filter; Iq_filter]
    #Update current
    current_r[1] += I_RI[1]
    current_i[1] += I_RI[2]
end
