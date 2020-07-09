function initialize_filter!(
    device_states,
    device::PSY.DynamicInverter{C, O, IC, DC, P, PSY.LCLFilter},
) where {
    C <: PSY.Converter,
    O <: PSY.OuterControl,
    IC <: PSY.InnerControl,
    DC <: PSY.DCSource,
    P <: PSY.FrequencyEstimator,
}
    #PowerFlow Data
    static_inj = PSY.get_static_injector(device)
    P0 = PSY.get_activepower(static_inj) / PSY.get_basepower(static_inj)
    Q0 = PSY.get_reactivepower(static_inj) / PSY.get_basepower(static_inj)
    Vm = PSY.get_voltage(PSY.get_bus(static_inj))
    θ = PSY.get_angle(PSY.get_bus(static_inj))
    S0 = P0 + Q0 * 1im
    V_R = Vm * cos(θ)
    V_I = Vm * sin(θ)
    V = V_R + V_I * 1im
    I = conj(S0 / V)
    I_R = real(I)
    I_I = imag(I)

    #Get Parameters
    filter = PSY.get_filter(device)
    lf = PSY.get_lf(filter)
    rf = PSY.get_rf(filter)
    cf = PSY.get_cf(filter)
    lg = PSY.get_lg(filter)
    rg = PSY.get_rg(filter)

    #Set parameters
    Ir_filter = I_R
    Ii_filter = I_I
    ω_sys = PSY.get_ω_ref(device)

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
        get_inner_vars(device)[VR_inv_var] = V_R
        get_inner_vars(device)[VI_inv_var] = V_I
        #Update Converter voltages
        get_inner_vars(device)[Vr_cnv_var] = sol_x0[1]
        get_inner_vars(device)[Vi_cnv_var] = sol_x0[2]
        #Update filter voltages
        get_inner_vars(device)[Vr_filter_var] = sol_x0[5]
        get_inner_vars(device)[Vi_filter_var] = sol_x0[6]
        #Update states
        filter_ix = get_local_state_ix(device, PSY.LCLFilter)
        filter_states = @view device_states[filter_ix]
        filter_states[1] = sol_x0[3] #Ir_cnv
        filter_states[2] = sol_x0[4] #Ii_cnv
        filter_states[3] = sol_x0[5] #Vr_filter
        filter_states[4] = sol_x0[6] #Vi_filter
        filter_states[5] = I_R #Ir_filter
        filter_states[6] = I_I #Ii_filter
    end
end
