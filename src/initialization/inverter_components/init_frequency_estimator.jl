function initialize_frequency_estimator!(
    device_states,
    static::PSY.StaticInjection,
    dyn_data::PSY.DynamicInverter{C, O, IC, DC, PSY.KauraPLL, F},
) where {
    C <: PSY.Converter,
    O <: PSY.OuterControl,
    IC <: PSY.InnerControl,
    DC <: PSY.DCSource,
    F <: PSY.Filter,
}
    Vr_filter = get_inner_vars(dyn_data)[Vr_filter_var]
    Vi_filter = get_inner_vars(dyn_data)[Vi_filter_var]

    θ0_pll = angle(Vr_filter + Vi_filter * 1im)
    Vpll_d0 = Vr_filter
    Vpll_q0 = 0.0
    ϵ_pll0 = 0.0

    pll_ix = get_local_state_ix(dyn_data, PSY.KauraPLL)

    #Obtain indices for component w/r to device
    local_ix = get_local_state_ix(dyn_data, PSY.KauraPLL)

    #Update guess of PLL states
    pll_states = @view device_states[local_ix]
    pll_states[1] = Vpll_d0
    pll_states[2] = Vpll_q0
    pll_states[3] = ϵ_pll0
    pll_states[4] = θ0_pll

    #Update guess of frequency estimator
    get_inner_vars(dyn_data)[ω_freq_estimator_var] = PSY.get_ω_ref(dyn_data)
end
