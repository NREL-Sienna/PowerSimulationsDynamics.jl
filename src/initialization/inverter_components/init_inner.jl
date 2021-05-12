function initialize_inner!(
    device_states,
    static::PSY.StaticInjection,
    dynamic_device::PSY.DynamicInverter{C, O, PSY.CurrentControl, DC, P, F},
) where {
    C <: PSY.Converter,
    O <: PSY.OuterControl,
    DC <: PSY.DCSource,
    P <: PSY.FrequencyEstimator,
    F <: PSY.Filter,
}

    #Obtain external states inputs for component
    external_ix = get_input_port_ix(dynamic_device, PSY.CurrentControl)
    Ir_filter = device_states[external_ix[1]]
    Ii_filter = device_states[external_ix[2]]
    Ir_cnv = device_states[external_ix[3]]
    Ii_cnv = device_states[external_ix[4]]
    Vr_filter = device_states[external_ix[5]] #TODO: Should be inner reference after initialization
    Vi_filter = device_states[external_ix[6]] #TODO: Should be inner reference after initialization

    #Obtain inner variables for component
    ω_oc = PSY.get_ω_ref(dynamic_device)
    θ0_oc = get_inner_vars(dynamic_device)[θ_oc_var]
    Vdc = get_inner_vars(dynamic_device)[Vdc_var]

    #Obtain output of converter
    Vr_cnv0 = get_inner_vars(dynamic_device)[Vr_cnv_var]
    Vi_cnv0 = get_inner_vars(dynamic_device)[Vi_cnv_var]

    #Get Voltage Controller parameters
    inner_control = PSY.get_inner_control(dynamic_device)
    filter = PSY.get_filter(dynamic_device)
    kpv = PSY.get_kpv(inner_control)
    kiv = PSY.get_kiv(inner_control)
    kffi = PSY.get_kffi(inner_control)
    cf = PSY.get_cf(filter)
    rv = PSY.get_rv(inner_control)
    lv = PSY.get_lv(inner_control)

    #Get Current Controller parameters
    kpc = PSY.get_kpc(inner_control)
    kic = PSY.get_kic(inner_control)
    kffv = PSY.get_kffv(inner_control)
    lf = PSY.get_lf(filter)
    ωad = PSY.get_ωad(inner_control)
    kad = PSY.get_kad(inner_control)

    function f!(out, x)
        θ_oc = x[1]
        v_refr = x[2]
        ξ_d = x[3]
        ξ_q = x[4]
        γ_d = x[5]
        γ_q = x[6]
        ϕ_d = x[7]
        ϕ_q = x[8]

        #Reference Frame Transformations
        I_dq_filter = ri_dq(θ_oc + pi / 2) * [Ir_filter; Ii_filter]
        I_dq_cnv = ri_dq(θ_oc + pi / 2) * [Ir_cnv; Ii_cnv]
        V_dq_filter = ri_dq(θ_oc + pi / 2) * [Vr_filter; Vi_filter]
        V_dq_cnv0 = ri_dq(θ_oc + pi / 2) * [Vr_cnv0; Vi_cnv0]

        #Voltage controller references
        Vd_filter_ref = (v_refr - rv * I_dq_filter[d] + ω_oc * lv * I_dq_filter[q])
        Vq_filter_ref = (-rv * I_dq_filter[q] - ω_oc * lv * I_dq_filter[d])

        #Current controller references
        Id_cnv_ref = (
            kpv * (Vd_filter_ref - V_dq_filter[d]) + kiv * ξ_d -
            cf * ω_oc * V_dq_filter[q] + kffi * I_dq_filter[d]
        )
        Iq_cnv_ref = (
            kpv * (Vq_filter_ref - V_dq_filter[q]) +
            kiv * ξ_q +
            cf * ω_oc * V_dq_filter[d] + kffi * I_dq_filter[q]
        )

        #References for Converter Output Voltage
        Vd_cnv_ref = (
            kpc * (Id_cnv_ref - I_dq_cnv[d]) + kic * γ_d - ω_oc * lf * I_dq_cnv[q] +
            kffv * V_dq_filter[d] - kad * (V_dq_filter[d] - ϕ_d)
        )
        Vq_cnv_ref = (
            kpc * (Iq_cnv_ref - I_dq_cnv[q]) + kic * γ_q + ω_oc * lf * I_dq_cnv[d] +
            kffv * V_dq_filter[q] - kad * (V_dq_filter[q] - ϕ_q)
        )

        out[1] = Vd_filter_ref - V_dq_filter[d]
        out[2] = Vq_filter_ref - V_dq_filter[q]
        out[3] = Id_cnv_ref - I_dq_cnv[d]
        out[4] = Iq_cnv_ref - I_dq_cnv[q]
        out[5] = V_dq_filter[d] - ϕ_d
        out[6] = V_dq_filter[q] - ϕ_q
        out[7] = Vd_cnv_ref - V_dq_cnv0[d]
        out[8] = Vq_cnv_ref - V_dq_cnv0[q]
    end
    x0 = [θ0_oc, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    sol = NLsolve.nlsolve(f!, x0)
    if !NLsolve.converged(sol)
        @warn("Initialization in Inner Control failed")
    else
        sol_x0 = sol.zero
        #Update angle:
        get_inner_vars(dynamic_device)[θ_oc_var] = sol_x0[1]
        outer_ix = get_local_state_ix(dynamic_device, O)
        outer_states = @view device_states[outer_ix]
        #Assumes that angle is in second position
        outer_states[2] = sol_x0[1]
        #Update V_ref (#TODO)
        PSY.get_ext(dynamic_device)[CONTROL_REFS][V_ref_index] = sol_x0[2]
        PSY.set_V_ref!(
            PSY.get_reactive_power(PSY.get_outer_control(dynamic_device)),
            sol_x0[2],
        )
        get_inner_vars(dynamic_device)[V_oc_var] = sol_x0[2]
        #Update Converter modulation
        m0_dq = (ri_dq(sol_x0[1] + pi / 2) * [Vr_cnv0; Vi_cnv0]) ./ Vdc
        get_inner_vars(dynamic_device)[md_var] = m0_dq[d]
        get_inner_vars(dynamic_device)[mq_var] = m0_dq[q]
        #Update states
        inner_ix = get_local_state_ix(dynamic_device, PSY.CurrentControl)
        inner_states = @view device_states[inner_ix]
        inner_states[1] = sol_x0[3] #ξ_d
        inner_states[2] = sol_x0[4] #ξ_q
        inner_states[3] = sol_x0[5] #γ_d
        inner_states[4] = sol_x0[6] #γ_q
        inner_states[5] = sol_x0[7] #ϕ_d
        inner_states[6] = sol_x0[8] #ϕ_q
    end
end
