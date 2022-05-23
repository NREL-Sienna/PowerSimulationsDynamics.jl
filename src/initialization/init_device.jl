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
    Zmag = R_th^2 + X_th^2

    function f!(out, x)
        V_R_internal = x[1]
        V_I_internal = x[2]

        out[1] =
            R_th * (V_R_internal - V_R) / Zmag + X_th * (V_I_internal - V_I) / Zmag - I_R
        out[2] =
            R_th * (V_I_internal - V_I) / Zmag - X_th * (V_R_internal - V_R) / Zmag - I_I
    end
    x0 = [V_R, V_I]
    sol = NLsolve.nlsolve(f!, x0)
    if !NLsolve.converged(sol)
        @warn("Initialization in Source failed")
    else
        sol_x0 = sol.zero
        #Update terminal voltages
        V_internal = sqrt(sol_x0[1]^2 + sol_x0[2]^2)
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
    Zmag = R_th^2 + X_th^2
    function f!(out, x)
        V_R_internal = x[1]
        V_I_internal = x[2]

        out[1] =
            R_th * (V_R_internal - V_R) / Zmag + X_th * (V_I_internal - V_I) / Zmag - I_R
        out[2] =
            R_th * (V_I_internal - V_I) / Zmag - X_th * (V_R_internal - V_R) / Zmag - I_I
    end
    x0 = [V_R, V_I]
    sol = NLsolve.nlsolve(f!, x0)
    if !NLsolve.converged(sol)
        @warn("Initialization in Periodic Variable Source failed")
    else
        sol_x0 = sol.zero
        #Update terminal voltages
        V_internal = sqrt(sol_x0[1]^2 + sol_x0[2]^2)
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
    I_R = R * (V_R_from - V_R_to) / Zmag + X * (V_I_from - V_I_to) / Zmag_sq
    I_I = R * (V_I_from - V_I_to) / Zmag - X * (V_R_from - V_R_to) / Zmag_sq
    #Update Current
    device_states[1] = I_R
    device_states[2] = I_I
    return device_states
end
