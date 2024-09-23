function mass_matrix_freq_estimator_entries!(
    mass_matrix,
    freq_estimator::P,
    global_index::Base.ImmutableDict{Symbol, Int64},
) where {P <: PSY.FrequencyEstimator}
    @debug "Using default mass matrix entries $P"
end

function mdl_freq_estimator_ode!(
    device_states::AbstractArray{<:ACCEPTED_REAL_TYPES},
    output_ode::AbstractArray{<:ACCEPTED_REAL_TYPES},
    p::AbstractArray{<:ACCEPTED_REAL_TYPES},
    inner_vars::AbstractArray{<:ACCEPTED_REAL_TYPES},
    ω_sys::ACCEPTED_REAL_TYPES,
    dynamic_device::DynamicWrapper{PSY.DynamicInverter{C, O, IC, DC, PSY.KauraPLL, F, L}},
    h,
    t,
) where {
    C <: PSY.Converter,
    O <: PSY.OuterControl,
    IC <: PSY.InnerControl,
    DC <: PSY.DCSource,
    F <: PSY.Filter,
    L <: Union{Nothing, PSY.OutputCurrentLimiter},
}

    #Obtain external states inputs for component
    external_ix = get_input_port_ix(dynamic_device, PSY.KauraPLL)
    Vr_filter = device_states[external_ix[1]]
    Vi_filter = device_states[external_ix[2]]

    #Get parameters
    params = p[:params][:FrequencyEstimator]
    ω_lp = params[:ω_lp]
    kp_pll = params[:kp_pll]
    ki_pll = params[:ki_pll]

    f0 = get_system_base_frequency(dynamic_device)
    ωb = 2.0 * pi * f0

    #Obtain indices for component w/r to device
    local_ix = get_local_state_ix(dynamic_device, PSY.KauraPLL)

    #Define internal states for frequency estimator
    internal_states = @view device_states[local_ix]
    vpll_d = internal_states[1]
    vpll_q = internal_states[2]
    ϵ_pll = internal_states[3]
    θ_pll = internal_states[4]

    #Transform to internal dq-PLL reference frame
    V_dq_pll = ri_dq(θ_pll + pi / 2) * [Vr_filter; Vi_filter]

    #Compute 6 states ODEs (D'Arco EPSR122 Model)
    #Output Voltage LPF (internal state)
    #𝜕vpll_d/𝜕t, D'Arco ESPR122 eqn. 12
    output_ode[local_ix[1]] = low_pass(V_dq_pll[d], vpll_d, 1.0, 1.0 / ω_lp)[2]
    #𝜕vpll_q/𝜕t, D'Arco ESPR122 eqn. 12
    output_ode[local_ix[2]] = low_pass(V_dq_pll[q], vpll_q, 1.0, 1.0 / ω_lp)[2]
    #PI Integrator (internal state)
    pi_output, dϵ_dt = pi_block(atan(vpll_q, vpll_d), ϵ_pll, kp_pll, ki_pll)
    #𝜕dϵ_pll/𝜕t, D'Arco ESPR122 eqn. 13
    output_ode[local_ix[3]] = dϵ_dt
    #PLL Frequency Deviation (internal state), Note: D'Arco ESPR122 eqn. 14 is missing (1.0-ω_sys) term. 
    #See Hug ISGT-EUROPE2018 Eqns. 26-28 for proper treatment of PLL reference frame. 
    Δω = 1.0 - ω_sys + pi_output
    #𝜕θ_pll/𝜕t, D'Arco ESPR122 eqn. 15
    output_ode[local_ix[4]] = ωb * Δω

    #Update inner_vars
    #PLL frequency, D'Arco EPSR122 eqn. 16
    inner_vars[ω_freq_estimator_var] = Δω + ω_sys
    inner_vars[θ_freq_estimator_var] = θ_pll
    return
end

function mdl_freq_estimator_ode!(
    device_states::AbstractArray{<:ACCEPTED_REAL_TYPES},
    output_ode::AbstractArray{<:ACCEPTED_REAL_TYPES},
    p::AbstractArray{<:ACCEPTED_REAL_TYPES},
    inner_vars::AbstractArray{<:ACCEPTED_REAL_TYPES},
    ω_sys::ACCEPTED_REAL_TYPES,
    dynamic_device::DynamicWrapper{
        PSY.DynamicInverter{C, O, IC, DC, PSY.ReducedOrderPLL, F, L},
    },
    h,
    t,
) where {
    C <: PSY.Converter,
    O <: PSY.OuterControl,
    IC <: PSY.InnerControl,
    DC <: PSY.DCSource,
    F <: PSY.Filter,
    L <: Union{Nothing, PSY.OutputCurrentLimiter},
}

    #Obtain external states inputs for component
    external_ix = get_input_port_ix(dynamic_device, PSY.ReducedOrderPLL)
    Vr_filter = device_states[external_ix[1]]
    Vi_filter = device_states[external_ix[2]]

    #Get parameters
    params = p[:params][:FrequencyEstimator]
    ω_lp = params[:ω_lp]
    kp_pll = params[:kp_pll]
    ki_pll = params[:ki_pll]
    f0 = get_system_base_frequency(dynamic_device)
    ωb = 2.0 * pi * f0

    #Obtain indices for component w/r to device
    local_ix = get_local_state_ix(dynamic_device, PSY.ReducedOrderPLL)

    #Define internal states for frequency estimator
    internal_states = @view device_states[local_ix]
    vpll_q = internal_states[1]
    ϵ_pll = internal_states[2]
    θ_pll = internal_states[3]

    #Transform to internal dq-PLL reference frame
    V_dq_pll = ri_dq(θ_pll + pi / 2) * [Vr_filter; Vi_filter]

    #Output Voltage LPF (internal state)
    #𝜕vpll_q/𝜕t, Low Pass Filter, Hug ISGT-EUROPE2018 eqn. 26
    output_ode[local_ix[1]] = low_pass(V_dq_pll[q], vpll_q, 1.0, 1.0 / ω_lp)[2]
    #PI Integrator (internal state)
    pi_output, dϵ_dt = pi_block(vpll_q, ϵ_pll, kp_pll, ki_pll)
    #𝜕dϵ_pll/𝜕t, Hug ISGT-EUROPE2018 eqn. 10
    output_ode[local_ix[2]] = dϵ_dt
    #PLL Frequency Deviation (internal state), Hug ISGT-EUROPE2018 eqn. 26 
    Δω = 1.0 - ω_sys + pi_output
    #𝜕θ_pll/𝜕t, Hug ISGT-EUROPE2018 eqns. 9, 26, 27 
    output_ode[local_ix[3]] = ωb * Δω

    #Update inner_vars
    #PLL frequency, D'Arco EPSR122 eqn. 16
    inner_vars[ω_freq_estimator_var] = Δω + ω_sys
    inner_vars[θ_freq_estimator_var] = θ_pll
    return
end

function mdl_freq_estimator_ode!(
    ::AbstractArray{<:ACCEPTED_REAL_TYPES},
    ::AbstractArray{<:ACCEPTED_REAL_TYPES},
    p::AbstractArray{<:ACCEPTED_REAL_TYPES},
    inner_vars::AbstractArray{<:ACCEPTED_REAL_TYPES},
    ω_sys::ACCEPTED_REAL_TYPES,
    dynamic_device::DynamicWrapper{
        PSY.DynamicInverter{C, O, IC, DC, PSY.FixedFrequency, F, L},
    },
    h,
    t,
) where {
    C <: PSY.Converter,
    O <: PSY.OuterControl,
    IC <: PSY.InnerControl,
    DC <: PSY.DCSource,
    F <: PSY.Filter,
    L <: Union{Nothing, PSY.OutputCurrentLimiter},
}
    #Get parameters
    frequency = p[:params][:FrequencyEstimator][:frequency]
    #Update inner_vars
    #PLL frequency
    inner_vars[ω_freq_estimator_var] = frequency
    return
end
