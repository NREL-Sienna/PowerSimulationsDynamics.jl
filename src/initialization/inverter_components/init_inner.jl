function initialize_inner!(
    device_states,
    p,
    static::PSY.StaticInjection,
    dynamic_device::DynamicWrapper{
        PSY.DynamicInverter{C, O, PSY.VoltageModeControl, DC, P, F, L},
    },
    inner_vars::AbstractVector,
) where {
    C <: PSY.Converter,
    O <: PSY.OuterControl,
    DC <: PSY.DCSource,
    P <: PSY.FrequencyEstimator,
    F <: PSY.Filter,
    L <: Union{Nothing, PSY.OutputCurrentLimiter},
}

    #Obtain external states inputs for component
    external_ix = get_input_port_ix(dynamic_device, PSY.VoltageModeControl)
    Ir_filter = device_states[external_ix[1]]
    Ii_filter = device_states[external_ix[2]]
    Ir_cnv = device_states[external_ix[3]]
    Ii_cnv = device_states[external_ix[4]]
    Vr_filter = device_states[external_ix[5]]
    Vi_filter = device_states[external_ix[6]]

    #Obtain inner variables for component
    ω_oc = p[:refs][:ω_ref]
    θ0_oc = inner_vars[θ_oc_var]
    Vdc = inner_vars[Vdc_var]

    #Obtain output of converter
    Vr_cnv0 = inner_vars[Vr_cnv_var]
    Vi_cnv0 = inner_vars[Vi_cnv_var]

    #Get Voltage Controller parameters
    params = p[:params][:InnerControl]
    cf = p[:params][:Filter][:cf]
    lf = p[:params][:Filter][:lf]
    function f!(out, x, params)
        θ_oc = x[1]
        v_refr = x[2]
        ξ_d = x[3]
        ξ_q = x[4]
        γ_d = x[5]
        γ_q = x[6]
        ϕ_d = x[7]
        ϕ_q = x[8]

        kpv = params[:kpv]
        kiv = params[:kiv]
        kffv = params[:kffv]
        rv = params[:rv]
        lv = params[:lv]
        kpc = params[:kpc]
        kic = params[:kic]
        kffi = params[:kffi]
        kad = params[:kad]

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
            cf * ω_oc * V_dq_filter[d] +
            kffi * I_dq_filter[q]
        )

        #References for Converter Output Voltage
        Vd_cnv_ref = (
            kpc * (Id_cnv_ref - I_dq_cnv[d]) + kic * γ_d - ω_oc * lf * I_dq_cnv[q] +
            kffv * V_dq_filter[d] - kad * (V_dq_filter[d] - ϕ_d)
        )
        Vq_cnv_ref = (
            kpc * (Iq_cnv_ref - I_dq_cnv[q]) +
            kic * γ_q +
            ω_oc * lf * I_dq_cnv[d] +
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
    prob = NonlinearSolve.NonlinearProblem{true}(f!, x0, params)
    sol = NonlinearSolve.solve(
        prob,
        NonlinearSolve.TrustRegion();
        sensealg = SciMLSensitivity.SteadyStateAdjoint(),
        reltol = STRICT_NLSOLVE_F_TOLERANCE,
        abstol = STRICT_NLSOLVE_F_TOLERANCE,
    )
    if !SciMLBase.successful_retcode(sol)
        @warn("Initialization in Inner Control failed")
    else
        sol_x0 = sol.u
        #Update angle:
        inner_vars[θ_oc_var] = sol_x0[1]
        outer_ix = get_local_state_ix(dynamic_device, O)
        outer_states = @view device_states[outer_ix]
        #Assumes that angle is in second position
        outer_states[1] = sol_x0[1]
        inner_vars[θ_oc_var] = sol_x0[1]
        set_V_ref!(p, sol_x0[2])
        PSY.set_V_ref!(
            PSY.get_reactive_power_control(PSY.get_outer_control(dynamic_device)),
            sol_x0[2],
        )
        inner_vars[V_oc_var] = sol_x0[2]
        # Update state if OuterControl is VOC
        if O <: PSY.OuterControl{PSY.ActiveVirtualOscillator, PSY.ReactiveVirtualOscillator}
            outer_states[2] = sol_x0[2]
        end
        #Update Converter modulation
        m0_dq = (ri_dq(sol_x0[1] + pi / 2) * [Vr_cnv0; Vi_cnv0]) ./ Vdc
        inner_vars[md_var] = m0_dq[d]
        inner_vars[mq_var] = m0_dq[q]
        #Update states
        inner_ix = get_local_state_ix(dynamic_device, PSY.VoltageModeControl)
        inner_states = @view device_states[inner_ix]
        inner_states[1] = sol_x0[3] #ξ_d
        inner_states[2] = sol_x0[4] #ξ_q
        inner_states[3] = sol_x0[5] #γ_d
        inner_states[4] = sol_x0[6] #γ_q
        inner_states[5] = sol_x0[7] #ϕ_d
        inner_states[6] = sol_x0[8] #ϕ_q
    end
    return
end

function initialize_inner!(
    device_states,
    p,
    static::PSY.StaticInjection,
    dynamic_device::DynamicWrapper{
        PSY.DynamicInverter{C, O, PSY.CurrentModeControl, DC, P, F, L},
    },
    inner_vars::AbstractVector,
) where {
    C <: PSY.Converter,
    O <: PSY.OuterControl,
    DC <: PSY.DCSource,
    P <: PSY.FrequencyEstimator,
    F <: PSY.Filter,
    L <: Union{Nothing, PSY.OutputCurrentLimiter},
}

    #Obtain external states inputs for component
    external_ix = get_input_port_ix(dynamic_device, PSY.CurrentModeControl)

    Ir_cnv = device_states[external_ix[3]]
    Ii_cnv = device_states[external_ix[4]]
    Vr_filter = device_states[external_ix[5]]
    Vi_filter = device_states[external_ix[6]]

    #Obtain inner variables for component
    ω_oc = p[:refs][:ω_ref]
    θ0_oc = inner_vars[θ_freq_estimator_var]
    Vdc = inner_vars[Vdc_var]
    Id_cnv_ref = inner_vars[Id_oc_var]
    Iq_cnv_ref = inner_vars[Iq_oc_var]

    #Obtain output of converter
    Vr_cnv0 = inner_vars[Vr_cnv_var]
    Vi_cnv0 = inner_vars[Vi_cnv_var]

    #Get Current Controller parameters
    kpc = p[:params][:InnerControl][:kpc]
    kic = p[:params][:InnerControl][:kic]
    kffv = p[:params][:InnerControl][:kffv]
    lf = p[:params][:Filter][:lf]
    params = [kpc, kic, kffv, lf]
    function f!(out, x, params)
        γ_d = x[1]
        γ_q = x[2]
        kpc, kic, kffv, lf = params
        #Reference Frame Transformations
        I_dq_cnv = ri_dq(θ0_oc + pi / 2) * [Ir_cnv; Ii_cnv]
        V_dq_filter = ri_dq(θ0_oc + pi / 2) * [Vr_filter; Vi_filter]
        V_dq_cnv0 = ri_dq(θ0_oc + pi / 2) * [Vr_cnv0; Vi_cnv0]

        #References for Converter Output Voltage
        Vd_cnv_ref = (
            kpc * (Id_cnv_ref - I_dq_cnv[d]) + kic * γ_d - ω_oc * lf * I_dq_cnv[q] +
            kffv * V_dq_filter[d]
        )
        Vq_cnv_ref = (
            kpc * (Iq_cnv_ref - I_dq_cnv[q]) +
            kic * γ_q +
            ω_oc * lf * I_dq_cnv[d] +
            kffv * V_dq_filter[q]
        )

        out[1] = Vd_cnv_ref - V_dq_cnv0[d]
        out[2] = Vq_cnv_ref - V_dq_cnv0[q]
    end
    x0 = [0.0, 0.0]
    prob = NonlinearSolve.NonlinearProblem{true}(f!, x0, params)
    sol = NonlinearSolve.solve(
        prob,
        NonlinearSolve.TrustRegion();
        sensealg = SciMLSensitivity.SteadyStateAdjoint(),
        reltol = STRICT_NLSOLVE_F_TOLERANCE,
        abstol = STRICT_NLSOLVE_F_TOLERANCE,
    )
    if !SciMLBase.successful_retcode(sol)
        @warn(
            "Initialization of AVR in $(PSY.get_name(static)) failed"
        )
    else
        sol_x0 = sol.u
        #Update Converter modulation
        m0_dq = (ri_dq(θ0_oc + pi / 2) * [Vr_cnv0; Vi_cnv0]) ./ Vdc
        inner_vars[md_var] = m0_dq[d]
        inner_vars[mq_var] = m0_dq[q]
        #Update states
        inner_ix = get_local_state_ix(dynamic_device, PSY.CurrentModeControl)
        inner_states = @view device_states[inner_ix]
        inner_states[1] = sol_x0[1] #γ_d
        inner_states[2] = sol_x0[2] #γ_q
    end
    return
end

function initialize_inner!(
    device_states,
    p,
    static::PSY.StaticInjection,
    dynamic_device::DynamicWrapper{
        PSY.DynamicInverter{C, O, PSY.RECurrentControlB, DC, P, F, L},
    },
    inner_vars::AbstractVector,
) where {
    C <: PSY.Converter,
    O <: PSY.OuterControl,
    DC <: PSY.DCSource,
    P <: PSY.FrequencyEstimator,
    F <: PSY.Filter,
    L <: Union{Nothing, PSY.OutputCurrentLimiter},
}
    # Obtain inner variables for component
    Vr_filter = inner_vars[Vr_filter_var]
    Vi_filter = inner_vars[Vi_filter_var]
    V_t = sqrt(Vr_filter^2 + Vi_filter^2)
    Iq_cmd = inner_vars[Iq_ic_var]
    Ip_oc = inner_vars[Id_oc_var]
    Iq_oc = inner_vars[Iq_oc_var]

    #Get Current Controller parameters
    inner_control = PSY.get_inner_control(dynamic_device)
    Q_Flag = PSY.get_Q_Flag(inner_control)
    PQ_Flag = PSY.get_PQ_Flag(inner_control)
    params = @view(@view(p[:params])[:InnerControl])
    V_ref0 = params[:V_ref0]
    K_vi = params[:K_vi]

    Ip_min, Ip_max, Iq_min, Iq_max =
        current_limit_logic(inner_control, Val(PQ_Flag), V_t, Ip_oc, Iq_cmd)

    if Ip_oc >= Ip_max + BOUNDS_TOLERANCE || Ip_min - BOUNDS_TOLERANCE >= Ip_oc
        @error(
            "Inverter $(PSY.get_name(static)) active current $(Ip_oc) out of limits $(Ip_min) $(Ip_max). Check Power Flow or Parameters"
        )
    end

    if Iq_oc >= Iq_max + BOUNDS_TOLERANCE || Iq_min - BOUNDS_TOLERANCE >= Iq_oc
        @error(
            "Inverter $(PSY.get_name(static)) reactive current $(Iq_oc) out of limits $(Iq_min) $(Iq_max). Check Power Flow or Parameters"
        )
    end

    if Q_Flag == 0
        local_ix = get_local_state_ix(dynamic_device, PSY.RECurrentControlB)
        #Define internal states for Inner Control
        internal_states = @view device_states[local_ix]
        internal_states[1] = V_t
        internal_states[2] = Iq_oc
    else
        local_ix = get_local_state_ix(dynamic_device, PSY.RECurrentControlB)
        K_vi = PSY.get_K_vi(inner_control)
        #Define internal states for Inner Control
        internal_states = @view device_states[local_ix]
        internal_states[1] = V_t
        internal_states[2] = Iq_cmd / K_vi
    end

    #Update additional variables
    # Based on PSS/E manual, if user does not provide V_ref0, then
    # V_ref0 is considered to be the output voltage of the PF solution
    if V_ref0 == 0.0
        params[:V_ref0] = V_t
    end
    return
end
