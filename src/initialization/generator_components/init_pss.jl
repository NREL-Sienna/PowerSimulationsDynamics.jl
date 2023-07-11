function initialize_pss!(
    device_states,
    static::PSY.StaticInjection,
    dynamic_device::DynamicWrapper{PSY.DynamicGenerator{M, S, A, TG, PSY.PSSFixed}},
    inner_vars::AbstractVector,
) where {M <: PSY.Machine, S <: PSY.Shaft, A <: PSY.AVR, TG <: PSY.TurbineGov} end

function initialize_pss!(
    device_states,
    static::PSY.StaticInjection,
    dynamic_device::DynamicWrapper{PSY.DynamicGenerator{M, S, A, TG, PSY.IEEEST}},
    inner_vars::AbstractVector,
) where {M <: PSY.Machine, S <: PSY.Shaft, A <: PSY.AVR, TG <: PSY.TurbineGov}
    #Get Signal Input Integer
    pss = PSY.get_pss(dynamic_device)
    input_code = PSY.get_input_code(pss)
    #Remote bus control not supported

    #Get Input Signal
    u = get_pss_input_signal(
        Val(input_code),
        device_states,
        inner_vars,
        1.0,
        dynamic_device,
    )

    #Obtain PSS States
    pss_ix = get_local_state_ix(dynamic_device, typeof(pss))
    pss_states = @view device_states[pss_ix]

    # Get Parameters
    A1 = PSY.get_A1(pss)
    A2 = PSY.get_A2(pss)
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

    #Error non-valid parameters
    if A6 > eps() && A2 < eps()
        error(
            "A6 cannot be greater than zero if A2 is zero. Correct data in PSS associated with $(PSY.get_name(dynamic_device))",
        )
    end

    if T1 > eps() && T2 < eps()
        error(
            "T2 cannot be greater than zero if T1 is zero. Correct data in PSS associated with $(PSY.get_name(dynamic_device))",
        )
    end

    if T3 > eps() && T4 < eps()
        error(
            "T3 cannot be greater than zero if T4 is zero. Correct data in PSS associated with $(PSY.get_name(dynamic_device))",
        )
    end

    if T5 < eps()
        error(
            "T5 is not allowed to be zero. Correct data in PSS associated with $(PSY.get_name(dynamic_device))",
        )
    end

    #Compute Parameter Ratios
    A6_A2 = A2 < eps() ? 0.0 : (A6 / A2)
    T1_T2 = T2 < eps() ? 0.0 : (T1 / T2)
    T3_T4 = T4 < eps() ? 0.0 : (T3 / T4)
    KsT5_T6 = T6 < eps() ? 0.0 : (Ks * T5 / T6)

    #Compute steady-state values
    x_p1 = 0.0
    x_p2 = u
    x_p3 = 0.0
    x_p4 = u
    y_f = A6_A2 * x_p2 + (A5 - A1 * A6_A2) * x_p3 + (1.0 - A6_A2) * x_p4
    x_p5 = y_f * (1.0 - T1_T2)
    y_LL1 = x_p5 + T1_T2 * y_f
    x_p6 = y_LL1 * (1.0 - T3_T4)
    y_LL2 = x_p6 + T3_T4 * y_LL1
    x_p7 = y_LL2
    y_out = KsT5_T6 * (y_LL2 - x_p7)

    #Compute and update output signal
    V_ss = clamp(y_out, Ls_min, Ls_max)
    #Compute compensated terminal voltage
    V_R = inner_vars[VR_gen_var]
    V_I = inner_vars[VI_gen_var]
    #To do: Figure out how to compensate terminal voltage
    V_ct = sqrt(V_R^2 + V_I^2)

    #Compute PSS output signal
    V_pss = output_pss_limiter(V_ss, V_ct, V_cl, V_cu)

    #Update Inner Vars
    inner_vars[V_pss_var] = V_pss

    #Update States
    pss_states[1] = x_p1
    pss_states[2] = x_p2
    pss_states[3] = x_p3
    pss_states[4] = x_p4
    pss_states[5] = x_p5
    pss_states[6] = x_p6
    pss_states[7] = x_p7
    return
end

function initialize_pss!(
    device_states,
    static::PSY.StaticInjection,
    dynamic_device::DynamicWrapper{PSY.DynamicGenerator{M, S, A, TG, PSY.STAB1}},
    inner_vars::AbstractVector,
) where {M <: PSY.Machine, S <: PSY.Shaft, A <: PSY.AVR, TG <: PSY.TurbineGov}
    #Get Signal Input Integer
    pss = PSY.get_pss(dynamic_device)

    #Obtain PSS States
    pss_ix = get_local_state_ix(dynamic_device, typeof(pss))
    pss_states = @view device_states[pss_ix]

    #Compute steady-state values
    x_p1 = 0.0
    x_p2 = 0.0
    x_p3 = 0.0

    #Compute PSS output signal
    V_pss = 0.0

    #Update Inner Vars
    inner_vars[V_pss_var] = V_pss

    #Update States
    pss_states[1] = x_p1
    pss_states[2] = x_p2
    pss_states[3] = x_p3
    return
end

function initialize_pss!(
    device_states,
    static::PSY.StaticInjection,
    dynamic_device::DynamicWrapper{PSY.DynamicGenerator{M, S, A, TG, PSY.PSS2A}},
    inner_vars::AbstractVector,
) where {M <: PSY.Machine, S <: PSY.Shaft, A <: PSY.AVR, TG <: PSY.TurbineGov}
    #Get Signal Input Integer
    pss = PSY.get_pss(dynamic_device)
    #Remote bus control not supported

    basepower = PSY.get_base_power(dynamic_device)
    Sbase = get_system_base_power(dynamic_device)

    #Obtain external states inputs for component
    external_ix = get_input_port_ix(dynamic_device, typeof(pss))
    ω = device_states[external_ix[1]]

    #Get Input Signal 1
    u_1 = get_pss_input_signal(
        Val(PSY.get_input_code_1(pss)),
        device_states,
        inner_vars,
        1.0,
        dynamic_device,
    )

    if PSY.get_input_code_1(pss) == 3
        u_1 = u_1 * (Sbase / basepower) * ω
    end

    #Get Input Signal 2
    u_2 = get_pss_input_signal(
        Val(PSY.get_input_code_2(pss)),
        device_states,
        inner_vars,
        1.0,
        dynamic_device,
    )

    if PSY.get_input_code_2(pss) == 3
        u_2 = u_2 * (Sbase / basepower) * ω
    end

    #Obtain PSS States
    pss_ix = get_local_state_ix(dynamic_device, typeof(pss))
    pss_states = @view device_states[pss_ix]

    # Get Required Parameters
    M_rtf = PSY.get_M_rtf(pss)
    N_rtf = PSY.get_N_rtf(pss)
    Tw1 = PSY.get_Tw1(pss)
    Tw3 = PSY.get_Tw3(pss)
    T9 = PSY.get_T9(pss)

    #Error non-valid parameters
    if M_rtf * N_rtf > 8
        error(
            "M*N cannot be greater than 8. Correct data in PSS associated with $(PSY.get_name(dynamic_device))",
        )
    end

    if Tw1 < eps()
        error(
            "Tw1 is not allowed to be zero. Correct data in PSS associated with $(PSY.get_name(dynamic_device))",
        )
    end

    if Tw3 < eps()
        error(
            "Tw3 is not allowed to be zero. Correct data in PSS associated with $(PSY.get_name(dynamic_device))",
        )
    end

    if T9 < eps()
        error(
            "T9 is not allowed to be zero. Correct data in PSS associated with $(PSY.get_name(dynamic_device))",
        )
    end

    #Compute steady-state values
    x_p1 = -u_1
    x_p2 = 0.0
    x_p3 = 0.0

    x_p4 = -u_2
    x_p5 = 0.0
    x_p6 = 0.0

    x_p7 = 0.0
    x_p8 = 0.0
    x_p9 = 0.0
    x_p10 = 0.0
    x_p11 = 0.0
    x_p12 = 0.0
    x_p13 = 0.0
    x_p14 = 0.0

    x_p15 = 0.0
    x_p16 = 0.0

    #Update Inner Vars
    inner_vars[V_pss_var] = 0.0

    #Update States
    pss_states[1] = x_p1
    pss_states[2] = x_p2
    pss_states[3] = x_p3
    pss_states[4] = x_p4
    pss_states[5] = x_p5
    pss_states[6] = x_p6
    pss_states[7] = x_p7
    pss_states[8] = x_p8
    pss_states[9] = x_p9
    pss_states[10] = x_p10
    pss_states[11] = x_p11
    pss_states[12] = x_p12
    pss_states[13] = x_p13
    pss_states[14] = x_p14
    pss_states[15] = x_p15
    pss_states[16] = x_p16
    return
end

function initialize_pss!(
    device_states,
    static::PSY.StaticInjection,
    dynamic_device::DynamicWrapper{PSY.DynamicGenerator{M, S, A, TG, PSY.PSS2B}},
    inner_vars::AbstractVector,
) where {M <: PSY.Machine, S <: PSY.Shaft, A <: PSY.AVR, TG <: PSY.TurbineGov}
    #Get Signal Input Integer
    pss = PSY.get_pss(dynamic_device)
    #Remote bus control not supported

    basepower = PSY.get_base_power(dynamic_device)
    Sbase = get_system_base_power(dynamic_device)

    #Obtain external states inputs for component
    external_ix = get_input_port_ix(dynamic_device, typeof(pss))
    ω = device_states[external_ix[1]]

    #Get Input Signal 1
    u_1 = get_pss_input_signal(
        Val(PSY.get_input_code_1(pss)),
        device_states,
        inner_vars,
        1.0,
        dynamic_device,
    )

    if PSY.get_input_code_1(pss) == 3
        u_1 = u_1 * (Sbase / basepower) * ω
    end

    #Get Input Signal 2
    u_2 = get_pss_input_signal(
        Val(PSY.get_input_code_2(pss)),
        device_states,
        inner_vars,
        1.0,
        dynamic_device,
    )

    if PSY.get_input_code_2(pss) == 3
        u_2 = u_2 * (Sbase / basepower) * ω
    end

    #Obtain PSS States
    pss_ix = get_local_state_ix(dynamic_device, typeof(pss))
    pss_states = @view device_states[pss_ix]

    # Get Required Parameters
    M_rtf = PSY.get_M_rtf(pss)
    N_rtf = PSY.get_N_rtf(pss)
    Tw1 = PSY.get_Tw1(pss)
    Tw3 = PSY.get_Tw3(pss)
    T9 = PSY.get_T9(pss)
    Vs1_min, Vs1_max = PSY.get_Vs1_lim(pss)
    Vs2_min, Vs2_max = PSY.get_Vs2_lim(pss)

    #Error for inputs outside limits
    if u_1 < Vs1_min || u_1 > Vs1_max
        error("First input for PSS $(PSY.get_name(dynamic_device)) is outside the limits")
    end

    if u_2 < Vs2_min || u_2 > Vs2_max
        error("Second input for PSS $(PSY.get_name(dynamic_device)) is outside the limits")
    end

    #Error non-valid parameters
    if M_rtf * N_rtf > 8
        error(
            "M*N cannot be greater than 8. Correct data in PSS associated with $(PSY.get_name(dynamic_device))",
        )
    end

    if Tw1 < eps()
        error(
            "Tw1 is not allowed to be zero. Correct data in PSS associated with $(PSY.get_name(dynamic_device))",
        )
    end

    if Tw3 < eps()
        error(
            "Tw3 is not allowed to be zero. Correct data in PSS associated with $(PSY.get_name(dynamic_device))",
        )
    end

    if T9 < eps()
        error(
            "T9 is not allowed to be zero. Correct data in PSS associated with $(PSY.get_name(dynamic_device))",
        )
    end

    #Compute steady-state values
    x_p1 = -u_1
    x_p2 = 0.0
    x_p3 = 0.0

    x_p4 = -u_2
    x_p5 = 0.0
    x_p6 = 0.0

    x_p7 = 0.0
    x_p8 = 0.0
    x_p9 = 0.0
    x_p10 = 0.0
    x_p11 = 0.0
    x_p12 = 0.0
    x_p13 = 0.0
    x_p14 = 0.0

    x_p15 = 0.0
    x_p16 = 0.0
    x_p17 = 0.0

    #Update Inner Vars
    inner_vars[V_pss_var] = 0.0

    #Update States
    pss_states[1] = x_p1
    pss_states[2] = x_p2
    pss_states[3] = x_p3
    pss_states[4] = x_p4
    pss_states[5] = x_p5
    pss_states[6] = x_p6
    pss_states[7] = x_p7
    pss_states[8] = x_p8
    pss_states[9] = x_p9
    pss_states[10] = x_p10
    pss_states[11] = x_p11
    pss_states[12] = x_p12
    pss_states[13] = x_p13
    pss_states[14] = x_p14
    pss_states[15] = x_p15
    pss_states[16] = x_p16
    pss_states[17] = x_p17
    return
end

#Currently not working properly.
#=
function initialize_pss!(
    device_states,
    static::PSY.StaticInjection,
    dynamic_device::DynamicWrapper{PSY.DynamicGenerator{M, S, A, TG, PSY.PSSSimple},
) where {M <: PSY.Machine, S <: PSY.Shaft, A <: PSY.AVR, TG <: PSY.TurbineGov}

    ω0 = get_ω_ref(dynamic_device)
    P0 = PSY.get_active_power(static)
    τe = get_inner_vars(dynamic_device)[τe_var]

    #Get parameters
    pss = PSY.get_pss(dynamic_device)
    K_ω = PSY.get_K_ω(pss)
    K_p = PSY.get_K_p(pss)

    #It assumes that the system is balanced so:
    get_inner_vars(dynamic_device)[V_pss_var] = K_ω * (ω0 - ω0) + K_p * (ω0 * τe - P0)

end
=#
