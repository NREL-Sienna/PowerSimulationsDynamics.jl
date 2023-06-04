function mass_matrix_pss_entries!(
    _,
    ::P,
    global_index::Base.ImmutableDict{Symbol, Int64},
) where {P <: PSY.PSS}
    @debug "Using default mass matrix entries $P"
end

function mass_matrix_pss_entries!(
    mass_matrix,
    pss::PSY.IEEEST,
    global_index::Base.ImmutableDict{Symbol, Int64},
)
    mass_matrix[global_index[:x_p1], global_index[:x_p1]] = PSY.get_A4(pss)
    mass_matrix[global_index[:x_p3], global_index[:x_p3]] = PSY.get_A2(pss)
    mass_matrix[global_index[:x_p5], global_index[:x_p5]] = PSY.get_T2(pss)
    mass_matrix[global_index[:x_p6], global_index[:x_p6]] = PSY.get_T4(pss)
    mass_matrix[global_index[:x_p7], global_index[:x_p7]] = PSY.get_T6(pss)
    return
end

#####################################################
############ Obtain PSS Input Signals  ##############
#####################################################

function get_pss_input_signal(
    ::Val{1},
    device_states::AbstractArray{<:ACCEPTED_REAL_TYPES},
    inner_vars::AbstractArray{<:ACCEPTED_REAL_TYPES},
    ω_sys::ACCEPTED_REAL_TYPES,
    dynamic_device::DynamicWrapper{PSY.DynamicGenerator{M, S, A, TG, P}},
) where {M <: PSY.Machine, S <: PSY.Shaft, A <: PSY.AVR, TG <: PSY.TurbineGov, P <: PSY.PSS}
    #Get PSS
    pss = PSY.get_pss(dynamic_device)
    #Obtain external states inputs for component
    external_ix = get_input_port_ix(dynamic_device, typeof(pss))
    ω = device_states[external_ix[1]]
    return ω - 1.0
end

function get_pss_input_signal(
    ::Val{2},
    device_states::AbstractArray{<:ACCEPTED_REAL_TYPES},
    inner_vars::AbstractArray{<:ACCEPTED_REAL_TYPES},
    ω_sys::ACCEPTED_REAL_TYPES,
    dynamic_device::DynamicWrapper{PSY.DynamicGenerator{M, S, A, TG, P}},
) where {M <: PSY.Machine, S <: PSY.Shaft, A <: PSY.AVR, TG <: PSY.TurbineGov, P <: PSY.PSS}
    # TODO: Frequency Input for PSS not properly supported yet"
    return ω_sys - 1.0
end

function get_pss_input_signal(
    ::Val{3},
    device_states::AbstractArray{<:ACCEPTED_REAL_TYPES},
    inner_vars::AbstractArray{<:ACCEPTED_REAL_TYPES},
    ω_sys::ACCEPTED_REAL_TYPES,
    dynamic_device::DynamicWrapper{PSY.DynamicGenerator{M, S, A, TG, P}},
) where {M <: PSY.Machine, S <: PSY.Shaft, A <: PSY.AVR, TG <: PSY.TurbineGov, P <: PSY.PSS}
    basepower = PSY.get_base_power(dynamic_device)
    Sbase = get_system_base_power(dynamic_device)
    return inner_vars[τe_var] * (basepower / Sbase)
end

function get_pss_input_signal(
    ::Val{4},
    device_states::AbstractArray{<:ACCEPTED_REAL_TYPES},
    inner_vars::AbstractArray{<:ACCEPTED_REAL_TYPES},
    ω_sys::ACCEPTED_REAL_TYPES,
    dynamic_device::DynamicWrapper{PSY.DynamicGenerator{M, S, A, TG, P}},
) where {M <: PSY.Machine, S <: PSY.Shaft, A <: PSY.AVR, TG <: PSY.TurbineGov, P <: PSY.PSS}
    P_ref = get_P_ref(device)
    return inner_vars[τm_var] - P_ref
end

function get_pss_input_signal(
    ::Val{5},
    device_states::AbstractArray{<:ACCEPTED_REAL_TYPES},
    inner_vars::AbstractArray{<:ACCEPTED_REAL_TYPES},
    ω_sys::ACCEPTED_REAL_TYPES,
    dynamic_device::DynamicWrapper{PSY.DynamicGenerator{M, S, A, TG, P}},
) where {M <: PSY.Machine, S <: PSY.Shaft, A <: PSY.AVR, TG <: PSY.TurbineGov, P <: PSY.PSS}
    V_R = inner_vars[VR_gen_var]
    V_I = inner_vars[VI_gen_var]
    return sqrt(V_R^2 + V_I^2)
end

function get_pss_input_signal(
    ::Val{6},
    device_states::AbstractArray{<:ACCEPTED_REAL_TYPES},
    inner_vars::AbstractArray{<:ACCEPTED_REAL_TYPES},
    ω_sys::ACCEPTED_REAL_TYPES,
    dynamic_device::DynamicWrapper{PSY.DynamicGenerator{M, S, A, TG, P}},
) where {M <: PSY.Machine, S <: PSY.Shaft, A <: PSY.AVR, TG <: PSY.TurbineGov, P <: PSY.PSS}
    error(
        "Derivative of voltage input in PSS at $(PSY.get_name(dynamic_device)) not supported",
    )
end

############################################
### ODE calculations via device dispatch ###
############################################

function mdl_pss_ode!(
    ::AbstractArray{<:ACCEPTED_REAL_TYPES},
    ::AbstractArray{<:ACCEPTED_REAL_TYPES},
    inner_vars::AbstractArray{<:ACCEPTED_REAL_TYPES},
    ω_sys::ACCEPTED_REAL_TYPES,
    dynamic_device::DynamicWrapper{PSY.DynamicGenerator{M, S, A, TG, PSY.PSSFixed}},
) where {M <: PSY.Machine, S <: PSY.Shaft, A <: PSY.AVR, TG <: PSY.TurbineGov}

    #Update V_pss on inner vars
    inner_vars[V_pss_var] = PSY.get_V_pss(PSY.get_pss(dynamic_device))

    return
end

function mdl_pss_ode!(
    device_states::AbstractArray{<:ACCEPTED_REAL_TYPES},
    output_ode::AbstractArray{<:ACCEPTED_REAL_TYPES},
    inner_vars::AbstractArray{<:ACCEPTED_REAL_TYPES},
    ω_sys::ACCEPTED_REAL_TYPES,
    dynamic_device::DynamicWrapper{PSY.DynamicGenerator{M, S, A, TG, PSY.IEEEST}},
) where {M <: PSY.Machine, S <: PSY.Shaft, A <: PSY.AVR, TG <: PSY.TurbineGov}

    #Get Signal Input Integer
    pss = PSY.get_pss(dynamic_device)
    #Remote bus control not supported

    #Get Input Signal
    u = get_pss_input_signal(
        Val(PSY.get_input_code(pss)),
        device_states,
        inner_vars,
        ω_sys,
        dynamic_device,
    )

    #Obtain indices for component w/r to device
    local_ix = get_local_state_ix(dynamic_device, PSY.IEEEST)

    #Define inner states for component
    internal_states = @view device_states[local_ix]
    x_p1 = internal_states[1]
    x_p2 = internal_states[2]
    x_p3 = internal_states[3]
    x_p4 = internal_states[4]
    x_p5 = internal_states[5]
    x_p6 = internal_states[6]
    x_p7 = internal_states[7]

    # Get Parameters
    A1 = PSY.get_A1(pss)
    A2 = PSY.get_A2(pss)
    A3 = PSY.get_A3(pss)
    A4 = PSY.get_A4(pss)
    A5 = PSY.get_A5(pss)
    A6 = PSY.get_A6(pss)
    T1 = PSY.get_T1(pss)
    T2 = PSY.get_T2(pss)
    T3 = PSY.get_T3(pss)
    T4 = PSY.get_T4(pss)
    T5 = PSY.get_T5(pss)
    T6 = PSY.get_T6(pss)
    Ks = PSY.get_Ks(pss)
    Ls_min, Ls_max = PSY.get_Ls_lim(pss)
    V_cu = PSY.get_Vcu(pss)
    V_cl = PSY.get_Vcl(pss)

    #Compute Parameter Ratios
    #A6_A2 = A2 < eps() ? 0.0 : (A6 / A2)
    #T1_T2 = T2 < eps() ? 0.0 : (T1 / T2)
    #T3_T4 = T4 < eps() ? 0.0 : (T3 / T4)
    #KsT5_T6 = T6 < eps() ? 0.0 : (Ks * T5 / T6)

    #Compute output of the filter
    #y_f = A6_A2 * x_p2 + (A5 - A1 * A6_A2) * x_p3 + (1.0 - A6_A2) * x_p4

    #Define Output of Lead-Lag blocks
    #y_LL1 = x_p5 + T1_T2 * y_f
    #y_LL2 = x_p6 + T3_T4 * y_LL1

    #Define Output of Feedback block
    #y_out = KsT5_T6 * (y_LL2 - x_p7)

    # Compute block derivatives
    _, dxp1_dt, dxp2_dt = low_pass_2nd_mass_matrix(u, x_p1, x_p2, 1.0, A3, A4)
    y_f, dxp3_dt, dxp4_dt = lead_lag_2nd_mass_matrix(x_p2, x_p3, x_p4, A1, A2, A5, A6)
    y_LL1, dxp5_dt = lead_lag_mass_matrix(y_f, x_p5, 1.0, T1, T2)
    y_LL2, dxp6_dt = lead_lag_mass_matrix(y_LL1, x_p6, 1.0, T3, T4)
    y_out, dxp7_dt = high_pass_mass_matrix(y_LL2, x_p7, Ks * T5, T6)

    #Compute 7 states PSS ODE
    output_ode[local_ix[1]] = dxp1_dt
    output_ode[local_ix[2]] = dxp2_dt
    output_ode[local_ix[3]] = dxp3_dt
    output_ode[local_ix[4]] = dxp4_dt
    output_ode[local_ix[5]] = dxp5_dt
    output_ode[local_ix[6]] = dxp6_dt
    output_ode[local_ix[7]] = dxp7_dt

    #Compute and update output signal
    V_ss = clamp(y_out, Ls_min, Ls_max)
    #Compute compensated terminal voltage
    V_R = inner_vars[VR_gen_var]
    V_I = inner_vars[VI_gen_var]
    #To do: Figure out how to compensate terminal voltage
    V_ct = sqrt(V_R^2 + V_I^2)

    #Compute PSS output signal and update inner vars
    inner_vars[V_pss_var] = output_pss_limiter(V_ss, V_ct, V_cl, V_cu)
    return
end

function mdl_pss_ode!(
    device_states::AbstractArray{<:ACCEPTED_REAL_TYPES},
    output_ode::AbstractArray{<:ACCEPTED_REAL_TYPES},
    inner_vars::AbstractArray{<:ACCEPTED_REAL_TYPES},
    ω_sys::ACCEPTED_REAL_TYPES,
    dynamic_device::DynamicWrapper{PSY.DynamicGenerator{M, S, A, TG, PSY.STAB1}},
) where {M <: PSY.Machine, S <: PSY.Shaft, A <: PSY.AVR, TG <: PSY.TurbineGov}

    #Get Signal Input Integer
    pss = PSY.get_pss(dynamic_device)

    # Get Input Signal - only speed deviation for STAB1
    external_ix = get_input_port_ix(dynamic_device, PSY.STAB1)
    ω = device_states[external_ix[1]] # get machine speed
    u = ω - 1.0

    #Obtain indices for component w/r to device
    local_ix = get_local_state_ix(dynamic_device, PSY.STAB1)

    #Define inner states for component
    internal_states = @view device_states[local_ix]
    x_p1 = internal_states[1] # state for washout filter
    x_p2 = internal_states[2] # state for Lead Lag 1
    x_p3 = internal_states[3] # state for Lead Lag 2

    # Get Parameters
    KT = PSY.get_KT(pss)
    T = PSY.get_T(pss)
    T1T3 = PSY.get_T1T3(pss)
    T3 = PSY.get_T3(pss)
    T2T4 = PSY.get_T2T4(pss)
    T4 = PSY.get_T4(pss)
    H_lim = PSY.get_H_lim(pss)

    K = T * KT
    T1 = T3 * T1T3
    T2 = T4 * T2T4

    # Compute block derivatives

    y_hp, dxp1_dt = high_pass(u, x_p1, K, T)
    y_LL1, dxp2_dt = lead_lag(y_hp, x_p2, 1.0, T1, T3)
    y_LL2, dxp3_dt = lead_lag(y_LL1, x_p3, 1.0, T2, T4)

    #Compute 3 states PSS ODE
    output_ode[local_ix[1]] = dxp1_dt
    output_ode[local_ix[2]] = dxp2_dt
    output_ode[local_ix[3]] = dxp3_dt

    #Compute and update output signal
    V_ss = clamp(y_LL2, -H_lim, H_lim)
    inner_vars[V_pss_var] = V_ss
    return
end

function mdl_pss_ode!(
    device_states::AbstractArray{<:ACCEPTED_REAL_TYPES},
    output_ode::AbstractArray{<:ACCEPTED_REAL_TYPES},
    inner_vars::AbstractArray{<:ACCEPTED_REAL_TYPES},
    ω_sys::ACCEPTED_REAL_TYPES,
    dynamic_device::DynamicWrapper{PSY.DynamicGenerator{M, S, A, TG, PSY.PSS2B}},
) where {M <: PSY.Machine, S <: PSY.Shaft, A <: PSY.AVR, TG <: PSY.TurbineGov}

    #Get Signal Input Integer
    pss = PSY.get_pss(dynamic_device)
    #Remote bus control not supported

    basepower = PSY.get_base_power(dynamic_device)
    Sbase = get_system_base_power(dynamic_device)

    #Get Input Signal 1
    u_1 = get_pss_input_signal(
        Val(PSY.get_input_code_1(pss)),
        device_states,
        inner_vars,
        ω_sys,
        dynamic_device,
    )

    if PSY.get_input_code_1(pss) == 3
        u_1 = u_1 * (Sbase / basepower)
    end

    #Get Input Signal 2
    u_2 = get_pss_input_signal(
        Val(PSY.get_input_code_2(pss)),
        device_states,
        inner_vars,
        ω_sys,
        dynamic_device,
    )

    if PSY.get_input_code_2(pss) == 3
        u_2 = u_2 * (Sbase / basepower)
    end

    #Obtain indices for component w/r to device
    local_ix = get_local_state_ix(dynamic_device, PSY.PSS2A)

    #Define inner states for component
    internal_states = @view device_states[local_ix]
    x_p1 = internal_states[1]
    x_p2 = internal_states[2]
    x_p3 = internal_states[3]
    x_p4 = internal_states[4]
    x_p5 = internal_states[5]
    x_p6 = internal_states[6]
    x_p7 = internal_states[7]
    x_p8 = internal_states[8]
    x_p9 = internal_states[9]
    x_p10 = internal_states[10]
    x_p11 = internal_states[11]
    x_p12 = internal_states[12]
    x_p13 = internal_states[13]
    x_p14 = internal_states[14]
    x_p15 = internal_states[15]
    x_p16 = internal_states[16]
    x_p17 = internal_states[17]

    # Get Parameters
    M_rtf = PSY.get_M_rtf(pss)
    N_rtf = PSY.get_N_rtf(pss)
    Tw1 = PSY.get_Tw1(pss)
    Tw2 = PSY.get_Tw2(pss)
    T6 = PSY.get_T6(pss)
    Tw3 = PSY.get_Tw3(pss)
    Tw4 = PSY.get_Tw4(pss)
    T7 = PSY.get_T7(pss)
    Ks2 = PSY.get_Ks2(pss)
    Ks3 = PSY.get_Ks3(pss)
    T8 = PSY.get_T8(pss)
    T9 = PSY.get_T9(pss)
    Ks1 = PSY.get_Ks1(pss)
    T1 = PSY.get_T1(pss)
    T2 = PSY.get_T2(pss)
    T3 = PSY.get_T3(pss)
    T4 = PSY.get_T4(pss)
    T10 = PSY.get_T10(pss)
    T11 = PSY.get_T11(pss)
    Vs1_min, Vs1_max = PSY.get_Vs1_lim(pss)
    Vs2_min, Vs2_max = PSY.get_Vs2_lim(pss)
    Vst_min, Vst_max = PSY.get_Vst_lim(pss)

    # Clamp inputs
    u1 = clamp(u1, Vs1_min, Vs1_max)
    u2 = clamp(u2, Vs2_min, Vs2_max)

    # Compute block derivatives
    y_w_1_1, dxp1_dt = high_pass(u_1, x_p1, Tw1, Tw1)
    y_w_2_1, dxp2_dt = high_pass(y_w_1_1, x_p2, Tw2, Tw2)

    # To bypass use T6 = 0.0
    yt_1 = y_w_2_1
    dxp3_dt = 0.0

    if T6 != 0.0
        yt_1, dxp3_dt = low_pass(y_w_2_1, x_p3, 1.0, T6)
    end

    y_w_1_2, dxp4_dt = high_pass(u_2, x_p4, Tw3, Tw3)
    y_w_2_2, dxp5_dt = high_pass(y_w_1_2, x_p5, Tw4, Tw4)

    # To bypass use T7 = 0.0
    yt_2 = y_w_2_2
    dxp6_dt = 0.0

    if T7 != 0.0
        yt_2, dxp6_dt = low_pass(y_w_2_2, x_p6, Ks2, T7)
    end

    # To bypass use M = N = 0
    y_rtf, dxp7_dt, dxp8_dt, dxp9_dt, dxp10_dt, dxp11_dt, dxp12_dt, dxp13_dt, dxp14_dt =
        ramp_tracking_filter(
            yt_1 + Ks3 * yt_2,
            x_p7,
            x_p8,
            x_p9,
            x_p10,
            x_p11,
            x_p12,
            x_p13,
            x_p14,
            T9,
            T8,
            M_rtf,
            N_rtf,
        )

    yll_1, dxp15_dt = lead_lag(Ks1 * (y_rtf - x_p6), x_p15, 1.0, T1, T2)
    yll_2, dxp16_dt = lead_lag(yll_1, x_p16, 1.0, T3, T4)
    yll_3, dxp17_dt = lead_lag(yll_2, x_p17, 1.0, T10, T11)

    #Compute 17 states PSS ODE
    output_ode[local_ix[1]] = dxp1_dt
    output_ode[local_ix[2]] = dxp2_dt
    output_ode[local_ix[3]] = dxp3_dt
    output_ode[local_ix[4]] = dxp4_dt
    output_ode[local_ix[5]] = dxp5_dt
    output_ode[local_ix[6]] = dxp6_dt
    output_ode[local_ix[7]] = dxp7_dt
    output_ode[local_ix[8]] = dxp8_dt
    output_ode[local_ix[9]] = dxp9_dt
    output_ode[local_ix[10]] = dxp10_dt
    output_ode[local_ix[11]] = dxp11_dt
    output_ode[local_ix[12]] = dxp12_dt
    output_ode[local_ix[13]] = dxp13_dt
    output_ode[local_ix[14]] = dxp14_dt
    output_ode[local_ix[15]] = dxp15_dt
    output_ode[local_ix[16]] = dxp16_dt
    output_ode[local_ix[17]] = dxp17_dt

    #Compute PSS output signal and update inner vars
    inner_vars[V_pss_var] = clamp(yll_3, Vst_min, Vst_max)
    return
end

#Currently not working properly.
#=
function mdl_pss_ode!(
    device_states,
    output_ode,
inner_vars,
    ω_sys,
    dynamic_device::DynamicWrapper{PSY.DynamicGenerator{M, S, A, TG, PSY.PSSSimple}},
) where {M <: PSY.Machine, S <: PSY.Shaft, A <: PSY.AVR, TG <: PSY.TurbineGov}

    #Get references
    P0 = PSY.get_active_power(static)

    #Obtain external states for device
    external_ix = get_input_port_ix(dynamic_device, PSY.PSSSimple)
    ω = @view device_states[external_ix]

    #Define external inner vars for component
    V_th = sqrt(get_inner_vars(dynamic_device)[VR_gen_var]^2 + get_inner_vars(dynamic_device)[VI_gen_var]^2)
    τe = get_inner_vars(dynamic_device)[τe_var]

    #Get parameters
    pss = PSY.get_pss(dynamic_device)
    K_ω = PSY.get_K_ω(pss)
    K_p = PSY.get_K_p(pss)

    #Update V_pss on inner vars
    get_inner_vars(dynamic_device)[V_pss_var] = K_ω * (ω[1] - ω_sys) + K_p * (ω[1] * τe - P0)

    return
end
=#
