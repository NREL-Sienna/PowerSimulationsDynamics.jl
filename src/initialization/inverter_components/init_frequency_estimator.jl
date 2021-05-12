function initialize_frequency_estimator!(
    device_states,
    static::PSY.StaticInjection,
    dynamic_device::PSY.DynamicInverter{C, O, IC, DC, PSY.KauraPLL, F},
) where {
    C <: PSY.Converter,
    O <: PSY.OuterControl,
    IC <: PSY.InnerControl,
    DC <: PSY.DCSource,
    F <: PSY.Filter,
}
    Vr_filter = get_inner_vars(dynamic_device)[Vr_filter_var]
    Vi_filter = get_inner_vars(dynamic_device)[Vi_filter_var]

    #Get parameters
    pll_control = PSY.get_freq_estimator(dynamic_device)
    ω_lp = PSY.get_ω_lp(pll_control)
    kp_pll = PSY.get_kp_pll(pll_control)
    ki_pll = PSY.get_ki_pll(pll_control)

    #Get initial guess
    θ0_pll = angle(Vr_filter + Vi_filter * 1im)
    Vpll_d0 = Vr_filter
    Vpll_q0 = 0.0
    ϵ_pll0 = 0.0

    function f!(out, x)
        vpll_d = x[1]
        vpll_q = x[2]
        ϵ_pll = x[3]
        θ_pll = x[4]

        V_dq_pll = ri_dq(θ_pll + pi / 2) * [Vr_filter; Vi_filter]

        out[1] = (V_dq_pll[d] - vpll_d)
        out[2] = (V_dq_pll[q] - vpll_q)
        out[3] = atan(vpll_q / vpll_d)
        out[4] = (kp_pll * atan(vpll_q / vpll_d) + ki_pll * ϵ_pll)
    end

    x0 = [Vpll_d0, Vpll_q0, ϵ_pll0, θ0_pll]
    sol = NLsolve.nlsolve(f!, x0)
    if !NLsolve.converged(sol)
        @warn("Initialization in PLL failed")
    else
        sol_x0 = sol.zero
        pll_ix = get_local_state_ix(dynamic_device, PSY.KauraPLL)

        #Obtain indices for component w/r to device
        local_ix = get_local_state_ix(dynamic_device, PSY.KauraPLL)

        #Update guess of PLL states
        pll_states = @view device_states[local_ix]
        pll_states[1] = sol_x0[1]
        pll_states[2] = sol_x0[2]
        pll_states[3] = sol_x0[3]
        pll_states[4] = sol_x0[4]

        #Update guess of frequency estimator
        get_inner_vars(dynamic_device)[ω_freq_estimator_var] = PSY.get_ω_ref(dynamic_device)
    end
end

function initialize_frequency_estimator!(
    device_states,
    static::PSY.StaticInjection,
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

    #Update guess of frequency estimator
    get_inner_vars(dynamic_device)[ω_freq_estimator_var] = frequency
end
