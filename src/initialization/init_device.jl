function initialize_dynamic_device!(
    dynamic_device::DynamicWrapper{DynG},
    static::PSY.StaticInjection,
    initial_inner_vars::AbstractVector,
) where {DynG <: PSY.DynamicGenerator}
    #Obtain States
    device_states = zeros(PSY.get_n_states(dynamic_device))

    #Initialize Machine and Shaft: δ and ω
    initialize_mach_shaft!(device_states, static, dynamic_device, initial_inner_vars)
    #Initialize extra Shaft states
    initialize_shaft!(device_states, static, dynamic_device, initial_inner_vars)
    #Initialize AVR
    initialize_avr!(device_states, static, dynamic_device, initial_inner_vars)
    #Initialize TG
    initialize_tg!(device_states, static, dynamic_device, initial_inner_vars)
    #Initialize PSS
    initialize_pss!(device_states, static, dynamic_device, initial_inner_vars)

    return device_states
end

function initialize_dynamic_device!(
    dynamic_device::DynamicWrapper{DynI},
    static::PSY.StaticInjection,
    initial_inner_vars::AbstractVector,
) where {DynI <: PSY.DynamicInverter}
    #Obtain States
    device_states = zeros(PSY.get_n_states(dynamic_device))

    #Initialize Machine and Shaft: V and I
    initialize_filter!(device_states, static, dynamic_device, initial_inner_vars)
    #Initialize Converter
    initialize_converter!(device_states, static, dynamic_device, initial_inner_vars)
    #Initialize freq estimator
    initialize_frequency_estimator!(
        device_states,
        static,
        dynamic_device,
        initial_inner_vars,
    )
    #Initialize OuterLoop
    initialize_outer!(device_states, static, dynamic_device, initial_inner_vars)
    #Initialize DCside
    initialize_DCside!(device_states, static, dynamic_device, initial_inner_vars)
    #Initialize Converter
    initialize_converter!(device_states, static, dynamic_device, initial_inner_vars)
    #Initialize InnerLoop
    initialize_inner!(device_states, static, dynamic_device, initial_inner_vars)
    return device_states
end

function initialize_static_device!(::StaticLoadWrapper)
    return
end

function initialize_static_device!(
    device::StaticWrapper{PSY.Source, T},
) where {T <: BusCategory}
    #PowerFlow Data
    P0 = PSY.get_active_power(device)
    Q0 = PSY.get_reactive_power(device)
    Vm = PSY.get_magnitude(PSY.get_bus(device))
    θ = PSY.get_angle(PSY.get_bus(device))
    S0 = P0 + Q0 * 1im
    V_R = Vm * cos(θ)
    V_I = Vm * sin(θ)
    V = V_R + V_I * 1im
    I = conj(S0 / V)
    I_R = real(I)
    I_I = imag(I)
    R_th = PSY.get_R_th(device.device)
    X_th = PSY.get_X_th(device.device)
    Zmag_sq_inv = 1 / (R_th^2 + X_th^2)

    function f!(out, x)
        V_R_internal = x[1]
        V_I_internal = x[2]

        out[1] =
            Zmag_sq_inv * (R_th * (V_R_internal - V_R) + X_th * (V_I_internal - V_I)) - I_R
        out[2] =
            Zmag_sq_inv * (R_th * (V_I_internal - V_I) - X_th * (V_R_internal - V_R)) - I_I
    end
    x0 = [V_R, V_I]
    sol = NLsolve.nlsolve(f!, x0)
    if !NLsolve.converged(sol)
        @warn("Initialization in Source failed")
    else
        sol_x0 = sol.zero
        #Update terminal voltages
        V_internal = hypot(sol_x0[1], sol_x0[2])
        θ_internal = atan(sol_x0[2], sol_x0[1])
        PSY.set_internal_voltage!(device.device, V_internal)
        PSY.set_internal_angle!(device.device, θ_internal)
        set_V_ref(device, PSY.get_internal_voltage(device.device))
        set_θ_ref(device, PSY.get_internal_angle(device.device))
    end
    return
end

function initialize_dynamic_device!(
    dynamic_device::DynamicWrapper{PSY.PeriodicVariableSource},
    source::PSY.Source,
    ::AbstractVector,
)
    device_states = zeros(PSY.get_n_states(dynamic_device))

    #PowerFlow Data
    P0 = PSY.get_active_power(source)
    Q0 = PSY.get_reactive_power(source)
    Vm = PSY.get_magnitude(PSY.get_bus(source))
    θ = PSY.get_angle(PSY.get_bus(source))
    S0 = P0 + Q0 * 1im
    V_R = Vm * cos(θ)
    V_I = Vm * sin(θ)
    V = V_R + V_I * 1im
    I = conj(S0 / V)
    I_R = real(I)
    I_I = imag(I)
    R_th = PSY.get_R_th(source)
    X_th = PSY.get_X_th(source)
    Zmag_sq_inv = 1 / (R_th^2 + X_th^2)
    function f!(out, x)
        V_R_internal = x[1]
        V_I_internal = x[2]

        out[1] =
            Zmag_sq_inv * (R_th * (V_R_internal - V_R) + X_th * (V_I_internal - V_I)) - I_R
        out[2] =
            Zmag_sq_inv * (R_th * (V_I_internal - V_I) - X_th * (V_R_internal - V_R)) - I_I
    end
    x0 = [V_R, V_I]
    sol = NLsolve.nlsolve(f!, x0)
    if !NLsolve.converged(sol)
        @warn("Initialization in Periodic Variable Source failed")
    else
        sol_x0 = sol.zero
        #Update terminal voltages
        V_internal = hypot(sol_x0[1], sol_x0[2])
        θ_internal = atan(sol_x0[2], sol_x0[1])

        V_internal_freqs = 0.0
        V_freqs = PSY.get_internal_voltage_frequencies(get_device(dynamic_device))
        V_coeff = PSY.get_internal_voltage_coefficients(get_device(dynamic_device))
        for (ix, ω) in enumerate(V_freqs)
            V_internal_freqs += V_coeff[ix][2]     #sin(0) = 0; cos(0)=1
        end

        θ_internal_freqs = 0.0
        θ_freqs = PSY.get_internal_angle_frequencies(get_device(dynamic_device))
        θ_coeff = PSY.get_internal_angle_coefficients(get_device(dynamic_device))
        for (ix, ω) in enumerate(θ_freqs)
            θ_internal_freqs += θ_coeff[ix][2]     #sin(0) = 0; cos(0)=1
        end
        V_internal_bias = V_internal - V_internal_freqs
        θ_internal_bias = θ_internal - θ_internal_freqs

        device_states[1] = V_internal
        device_states[2] = θ_internal
        PSY.set_internal_voltage_bias!(get_device(dynamic_device), V_internal_bias)
        PSY.set_internal_angle_bias!(get_device(dynamic_device), θ_internal_bias)
    end
    return device_states
end

function initialize_dynamic_device!(branch::BranchWrapper)
    device_states = zeros(PSY.get_n_states(branch))
    #PowerFlow Data
    arc = PSY.get_arc(branch)
    Vm_from = PSY.get_magnitude(PSY.get_from(arc))
    θ_from = PSY.get_angle(PSY.get_from(arc))
    Vm_to = PSY.get_magnitude(PSY.get_to(arc))
    θ_to = PSY.get_angle(PSY.get_to(arc))
    V_R_from = Vm_from * cos(θ_from)
    V_I_from = Vm_from * sin(θ_from)
    V_R_to = Vm_to * cos(θ_to)
    V_I_to = Vm_to * sin(θ_to)
    R = PSY.get_r(branch)
    X = PSY.get_x(branch)
    Zmag_sq = R^2 + X^2
    #Compute Current
    I_R = R * (V_R_from - V_R_to) / Zmag_sq + X * (V_I_from - V_I_to) / Zmag_sq
    I_I = R * (V_I_from - V_I_to) / Zmag_sq - X * (V_R_from - V_R_to) / Zmag_sq
    #Update Current
    device_states[1] = I_R
    device_states[2] = I_I
    return device_states
end

function initialize_dynamic_device!(
    dynamic_wrapper::DynamicWrapper{PSY.SingleCageInductionMachine},
    device::PSY.StaticInjection,
    ::AbstractVector,
)
    Sbase = get_system_base_power(dynamic_wrapper)

    #Obtain States
    device_states = zeros(PSY.get_n_states(dynamic_wrapper))

    # Get parameters
    dynamic_device = get_device(dynamic_wrapper)
    R_s = PSY.get_Rs(dynamic_device)
    X_ls = PSY.get_Xs(dynamic_device)
    X_m = PSY.get_Xm(dynamic_device)
    R_r = PSY.get_Rr(dynamic_device)
    X_lr = PSY.get_Xr(dynamic_device)
    A = PSY.get_A(dynamic_device)
    B = PSY.get_B(dynamic_device)
    C = PSY.get_C(dynamic_device)
    base_power = PSY.get_base_power(dynamic_device)

    #PowerFlow Data
    P0 = PSY.get_active_power(device) * Sbase / base_power # in pu (motor base)
    Q0 = PSY.get_reactive_power(device) * Sbase / base_power # in pu (motor base)
    Vm = PSY.get_magnitude(PSY.get_bus(device))
    θ = PSY.get_angle(PSY.get_bus(device))
    S0 = P0 + Q0 * 1im
    V_R = Vm * cos(θ)
    V_I = Vm * sin(θ)
    V = V_R + V_I * 1im
    I = conj(S0 / V) # total current (includes B_shunt + motor)
    I_R = real(I)
    I_I = imag(I)

    # voltages in qd reference frame
    v_ds = V_R
    v_qs = V_I

    # Reactances for the model
    X_ad = (1.0 / X_m + 1.0 / X_ls + 1.0 / X_lr)^(-1)
    X_aq = X_ad

    # Initial guess for NLSolve (assume B_shunt = 0)
    i_qs0 = I_I
    i_ds0 = I_R
    B_sh0 = 0.0
    ψ_qs0 = -v_ds + R_s * i_ds0
    ψ_ds0 = v_qs - R_s * i_qs0
    ψ_mq0 = ψ_qs0 - i_qs0 * X_ls
    ψ_md0 = ψ_ds0 - i_ds0 * X_ls
    ψ_qr0 = (ψ_mq0 / X_ad - ψ_qs0 / X_ls) * X_lr
    ψ_dr0 = (ψ_md0 / X_ad - ψ_ds0 / X_ls) * X_lr
    ωr0 = 0.98 # good guess for ind. motor
    τ_m00 = P0 / ωr0
    x0 = [i_qs0, i_ds0, B_sh0, ψ_qs0, ψ_ds0, ψ_qr0, ψ_dr0, ωr0, τ_m00]

    #
    function f!(out, x)
        i_qs = x[1]
        i_ds = x[2]
        B_sh = x[3]
        ψ_qs = x[4]
        ψ_ds = x[5]
        ψ_qr = x[6]
        ψ_dr = x[7]
        ωr = x[8]
        τ_m0 = x[9]
        ψ_mq = ψ_qs - i_qs * X_ls
        ψ_md = ψ_ds - i_ds * X_ls
        out[1] = -I_R + i_ds - V_I * B_sh # network interface
        out[2] = -I_I + i_qs + V_R * B_sh # network interface
        out[3] = v_qs - ψ_ds - R_s * i_qs #
        out[4] = v_ds + ψ_qs - R_s * i_ds #
        out[5] = -ψ_mq + X_aq * (ψ_qs / X_ls + ψ_qr / X_lr) # dψ_qs/dt = 0
        out[6] = -ψ_md + X_ad * (ψ_ds / X_ls + ψ_dr / X_lr) # dψ_ds/dt = 0
        out[7] = -(1 - ωr) * ψ_dr + R_r / X_lr * (ψ_mq - ψ_qr) # dψ_qr/dt = 0
        out[8] = (1 - ωr) * ψ_qr + R_r / X_lr * (ψ_md - ψ_dr) # dψ_dr/dt = 0
        out[9] = ψ_ds * i_qs - ψ_qs * i_ds - τ_m0 * (A * ωr^2 + B * ωr + C) # dωr/dt = 0
    end
    sol = NLsolve.nlsolve(f!, x0, ftol = STRICT_NLSOLVE_F_TOLERANCE)
    if !NLsolve.converged(sol)
        @warn("Initialization in Ind. Motor $(PSY.get_name(device)) failed")
    else
        sol_x0 = sol.zero
        device_states[1] = sol_x0[4] # ψ_qs
        device_states[2] = sol_x0[5] # ψ_ds
        device_states[3] = sol_x0[6] # ψ_qr
        device_states[4] = sol_x0[7] # ψ_dr
        device_states[5] = sol_x0[8] # ωr
        # update τ_ref and B_sh
        PSY.set_B_shunt!(dynamic_device, sol_x0[3]) # B_sh
        #set_B_shunt(dynamic_device, sol_x0[3]) # B_sh
        PSY.set_τ_ref!(dynamic_device, sol_x0[9]) # τ_m0
        #set_τ_ref(dynamic_device, sol_x0[9]) # τ_m0
    end
    return device_states
end
