function initialize_filter!(
    device_states,
    static::PSY.StaticInjection,
    dynamic_device::DynamicWrapper{PSY.DynamicInverter{C, O, IC, DC, P, PSY.LCLFilter}},
    inner_vars::AbstractVector,
) where {
    C <: PSY.Converter,
    O <: PSY.OuterControl,
    IC <: PSY.InnerControl,
    DC <: PSY.DCSource,
    P <: PSY.FrequencyEstimator,
}
    #PowerFlow Data
    P0 = PSY.get_active_power(static)
    Q0 = PSY.get_reactive_power(static)
    Vm = PSY.get_magnitude(PSY.get_bus(static))
    θ = PSY.get_angle(PSY.get_bus(static))
    S0 = P0 + Q0 * 1im
    V_R = Vm * cos(θ)
    V_I = Vm * sin(θ)
    V = V_R + V_I * 1im
    I = conj(S0 / V)
    I_R = real(I)
    I_I = imag(I)

    #Get Parameters
    filter = PSY.get_filter(dynamic_device)
    lf = PSY.get_lf(filter)
    rf = PSY.get_rf(filter)
    cf = PSY.get_cf(filter)
    lg = PSY.get_lg(filter)
    rg = PSY.get_rg(filter)

    #Set parameters
    Ir_filter = I_R
    Ii_filter = I_I
    ω_sys = get_ω_ref(dynamic_device)

    #To solve Vr_cnv, Vi_cnv, Ir_cnv, Ii_cnv, Vr_filter, Vi_filter
    function f!(out, x)
        Vr_cnv = x[1]
        Vi_cnv = x[2]
        Ir_cnv = x[3]
        Ii_cnv = x[4]
        Vr_filter = x[5]
        Vi_filter = x[6]

        #𝜕Ir_cnv/𝜕t
        out[1] = Vr_cnv - Vr_filter - rf * Ir_cnv + ω_sys * lf * Ii_cnv
        #𝜕Ii_cnv/𝜕t
        out[2] = Vi_cnv - Vi_filter - rf * Ii_cnv - ω_sys * lf * Ir_cnv
        #𝜕Vr_filter/𝜕t
        out[3] = Ir_cnv - Ir_filter + ω_sys * cf * Vi_filter
        #𝜕Vi_filter/𝜕t
        out[4] = Ii_cnv - Ii_filter - ω_sys * cf * Vr_filter
        #𝜕Ir_filter/𝜕t
        out[5] = Vr_filter - V_R - rg * Ir_filter + ω_sys * lg * Ii_filter
        #𝜕Ii_filter/𝜕t
        out[6] = Vi_filter - V_I - rg * Ii_filter - ω_sys * lg * Ir_filter
    end
    x0 = [V_R, V_I, I_R, I_I, V_R, V_I]
    sol = NLsolve.nlsolve(f!, x0)
    if !NLsolve.converged(sol)
        @warn("Initialization in Filter failed")
    else
        sol_x0 = sol.zero
        #Update terminal voltages
        inner_vars[Vr_inv_var] = V_R
        inner_vars[Vi_inv_var] = V_I
        #Update Converter voltages
        inner_vars[Vr_cnv_var] = sol_x0[1]
        inner_vars[Vi_cnv_var] = sol_x0[2]
        #Update filter voltages
        inner_vars[Vr_filter_var] = sol_x0[5]
        inner_vars[Vi_filter_var] = sol_x0[6]
        #Update states
        filter_ix = get_local_state_ix(dynamic_device, PSY.LCLFilter)
        filter_states = @view device_states[filter_ix]
        filter_states[1] = sol_x0[3] #Ir_cnv
        filter_states[2] = sol_x0[4] #Ii_cnv
        filter_states[3] = sol_x0[5] #Vr_filter
        filter_states[4] = sol_x0[6] #Vi_filter
        filter_states[5] = I_R #Ir_filter
        filter_states[6] = I_I #Ii_filter
    end
    return
end

function initialize_filter!(
    device_states,
    static::PSY.StaticInjection,
    dynamic_device::DynamicWrapper{PSY.DynamicInverter{C, O, IC, DC, P, PSY.RLFilter}},
    inner_vars::AbstractVector,
) where {
    C <: PSY.Converter,
    O <: PSY.OuterControl,
    IC <: PSY.InnerControl,
    DC <: PSY.DCSource,
    P <: PSY.FrequencyEstimator,
}
    #PowerFlow Data
    P0 = PSY.get_active_power(static)
    Q0 = PSY.get_reactive_power(static)
    Vm = PSY.get_magnitude(PSY.get_bus(static))
    θ = PSY.get_angle(PSY.get_bus(static))
    S0 = P0 + Q0 * 1im
    V_R = Vm * cos(θ)
    V_I = Vm * sin(θ)
    V = V_R + V_I * 1im
    # Doesn't match PSS/e initialization because PTI doesn't use the conjugate.
    I = conj(S0 / V)
    # PSS/e names I_R as Ip. But is calculated as P/Vt
    I_R = real(I)
    # PSS/e names I_I as Iq. But is calculated as Q/Vt
    I_I = imag(I)

    #Get Parameters
    filt = PSY.get_filter(dynamic_device)
    rf = PSY.get_rf(filt)
    lf = PSY.get_lf(filt)
    converter = PSY.get_converter(dynamic_device)
    R_source = PSY.get_R_source(converter)
    X_source = PSY.get_X_source(converter)

    #Update terminal voltages
    inner_vars[Vr_inv_var] = V_R
    inner_vars[Vi_inv_var] = V_I
    #Update filter currents (output of converter)
    inner_vars[Ir_inv_var] = I_R
    inner_vars[Ii_inv_var] = I_I

    #Update converter currents
    V_cnv = V + (rf + lf * 1im) * I
    I_aux = V_cnv / (R_source + X_source * 1im)
    I_cnv = I + I_aux

    # Update Control References
    PSY.set_Q_ref!(PSY.get_converter(dynamic_device), Q_bad)
    set_Q_ref(dynamic_device, Q_bad)
    PSY.set_Q_ref!(PSY.get_reactive_power(PSY.get_outer_control(dynamic_device)), Q_bad)
    PSY.set_P_ref!(PSY.get_active_power(PSY.get_outer_control(dynamic_device)), P0)
    set_P_ref(dynamic_device, P0)
    PSY.set_V_ref!(PSY.get_reactive_power(PSY.get_outer_control(dynamic_device)), Vm)
    set_V_ref(dynamic_device, Vm)

    inner_vars[Vr_cnv_var] = real(V_cnv)
    inner_vars[Vi_cnv_var] = imag(V_cnv)
    inner_vars[Ir_cnv_var] = real(I_cnv)
    inner_vars[Ii_cnv_var] = imag(I_cnv)
    return
end
