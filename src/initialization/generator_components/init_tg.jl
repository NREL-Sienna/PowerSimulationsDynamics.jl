function initialize_tg!(
    device_states,
    p,
    ::PSY.StaticInjection,
    dynamic_device::DynamicWrapper{PSY.DynamicGenerator{M, S, A, PSY.TGFixed, P}},
    inner_vars::AbstractVector,
) where {M <: PSY.Machine, S <: PSY.Shaft, A <: PSY.AVR, P <: PSY.PSS}
    tg = PSY.get_prime_mover(dynamic_device)
    eff = p[:params][:TurbineGov][:efficiency]
    τm0 = inner_vars[τm_var]

    P_ref = τm0 / eff
    #Update Control Refs
    PSY.set_P_ref!(tg, P_ref)
    set_P_ref!(p, P_ref)
    return
end

function initialize_tg!(
    device_states,
    device_parameters,
    static::PSY.StaticInjection,
    dynamic_device::DynamicWrapper{PSY.DynamicGenerator{M, S, A, PSY.TGTypeI, P}},
    inner_vars::AbstractVector,
) where {M <: PSY.Machine, S <: PSY.Shaft, A <: PSY.AVR, P <: PSY.PSS}

    #Get mechanical torque to SyncMach
    τm0 = inner_vars[τm_var]

    tg = PSY.get_prime_mover(dynamic_device)
    #Get Parameters
    local_ix_params = get_local_parameter_ix(dynamic_device, PSY.TGTypeI)
    internal_params = @view device_parameters[local_ix_params]
    R, _, Tc, T3, T4, T5, V_min, V_max = internal_params
    inv_R = R < eps() ? 0.0 : (1.0 / R)

    #Get References
    ω_ref = get_ω_ref(dynamic_device)
    ω0 = 1.0
    function f!(out, x)
        P_ref = x[1]
        x_g1 = x[2]
        x_g2 = x[3]
        x_g3 = x[4]

        #Compute auxiliary parameters
        P_in = P_ref + inv_R * (ω_ref - ω0)

        out[1] = P_in - x_g1
        out[2] = (1.0 - T3 / Tc) * x_g1 - x_g2
        out[3] = (1.0 - T4 / T5) * (x_g2 + (T3 / Tc) * x_g1) - x_g3
        out[4] = x_g3 + (T4 / T5) * (x_g2 + (T3 / Tc) * x_g1) - τm0
    end
    x0 = [
        τm0,
        τm0,
        (1.0 - T3 / Tc) * τm0,
        (1.0 - T4 / T5) * ((1.0 - T3 / Tc) * τm0 + (T3 / Tc) * τm0),
    ]
    sol = NLsolve.nlsolve(f!, x0; ftol = STRICT_NLSOLVE_F_TOLERANCE)
    if !NLsolve.converged(sol)
        CRC.@ignore_derivatives @warn(
            "Initialization of Turbine Governor $(PSY.get_name(static)) failed"
        )
    else
        sol_x0 = sol.zero
        #Update Control Refs
        PSY.set_P_ref!(tg, sol_x0[1])
        set_P_ref(dynamic_device, sol_x0[1])
        #Update states
        tg_ix = get_local_state_ix(dynamic_device, PSY.TGTypeI)
        tg_states = @view device_states[tg_ix]
        if (sol_x0[2] > V_max) || (sol_x0[2] < V_min)
            CRC.@ignore_derivatives @error(
                "Valve limits for TG in $(PSY.get_name(dynamic_device)) (x_g1 = $(sol_x0[2])), outside its limits V_max = $V_max, Vmin = $V_min. Consider updating the operating point."
            )
        end
        tg_states[1] = sol_x0[2]
        tg_states[2] = sol_x0[3]
        tg_states[3] = sol_x0[4]
    end
    return
end

function initialize_tg!(
    device_states,
    device_parameters,
    static::PSY.StaticInjection,
    dynamic_device::DynamicWrapper{PSY.DynamicGenerator{M, S, A, PSY.TGTypeII, P}},
    inner_vars::AbstractVector,
) where {M <: PSY.Machine, S <: PSY.Shaft, A <: PSY.AVR, P <: PSY.PSS}

    #Get mechanical torque to SyncMach
    τm0 = inner_vars[τm_var]
    #Get parameters
    tg = PSY.get_prime_mover(dynamic_device)
    local_ix_params = get_local_parameter_ix(dynamic_device, PSY.TGTypeII)
    internal_params = @view device_parameters[local_ix_params]
    R, T1, T2 = internal_params
    inv_R = R < eps() ? 0.0 : (1.0 / R)
    ω_ref = get_ω_ref(dynamic_device)
    ω0 = ω_ref

    function f!(out, x)
        P_ref = x[1]
        xg = x[2]
        out[1] = inv_R * (T1 / T2) * (ω_ref - ω0) + P_ref / 1.0 + xg - τm0
        out[2] = (1.0 / T2) * (inv_R * (1 - T2 / T1) * (ω_ref - ω0) - xg)
    end
    x0 = [τm0, 0.0]
    sol = NLsolve.nlsolve(f!, x0; ftol = STRICT_NLSOLVE_F_TOLERANCE)
    if !NLsolve.converged(sol)
        CRC.@ignore_derivatives @warn(
            "Initialization of Turbine Governor $(PSY.get_name(static)) failed"
        )
    else
        sol_x0 = sol.zero
        #Update Control Refs
        PSY.set_P_ref!(tg, sol_x0[1])
        set_P_ref(dynamic_device, sol_x0[1])
        #Update states
        tg_ix = get_local_state_ix(dynamic_device, PSY.TGTypeII)
        tg_states = @view device_states[tg_ix]
        tg_states[1] = sol_x0[2]
    end
    return
end

function initialize_tg!(
    device_states,
    device_parameters,
    static::PSY.StaticInjection,
    dynamic_device::DynamicWrapper{PSY.DynamicGenerator{M, S, A, PSY.GasTG, P}},
    inner_vars::AbstractVector,
) where {M <: PSY.Machine, S <: PSY.Shaft, A <: PSY.AVR, P <: PSY.PSS}

    #Get mechanical torque to SyncMach
    τm0 = inner_vars[τm_var]
    Δω = 0.0

    tg = PSY.get_prime_mover(dynamic_device)
    #Get parameters
    local_ix_params = get_local_parameter_ix(dynamic_device, PSY.GasTG)
    internal_params = @view device_parameters[local_ix_params]
    R, _, _, _, AT, Kt, V_min, V_max, D_turb = internal_params
    inv_R = R < eps() ? 0.0 : (1.0 / R)

    function f!(out, x)
        P_ref = x[1]
        x_g1 = x[2]
        x_g2 = x[3]
        x_g3 = x[4]

        x_in = min((P_ref - inv_R * Δω), (AT + Kt * (AT - x_g3)))
        out[1] = (x_in - x_g1) #dx_g1/dt
        x_g1_sat = V_min < x_g1 < V_max ? x_g1 : max(V_min, min(V_max, x_g1))
        out[2] = (x_g1_sat - x_g2)
        out[3] = (x_g2 - x_g3)
        out[4] = (x_g2 - D_turb * Δω) - τm0
    end
    x0 = [1.0 / inv_R, τm0, τm0, τm0]
    sol = NLsolve.nlsolve(f!, x0; ftol = STRICT_NLSOLVE_F_TOLERANCE)
    if !NLsolve.converged(sol)
        CRC.@ignore_derivatives @error(
            "Initialization of Turbine Governor $(PSY.get_name(static)) failed"
        )
    else
        sol_x0 = sol.zero
        #Update Control Refs
        PSY.set_P_ref!(tg, sol_x0[1])
        set_P_ref(dynamic_device, sol_x0[1])
        #Update states
        tg_ix = get_local_state_ix(dynamic_device, typeof(tg))
        tg_states = @view device_states[tg_ix]
        tg_states[1] = sol_x0[2]
        tg_states[2] = sol_x0[3]
        tg_states[3] = sol_x0[4]
    end
    return
end

function initialize_tg!(
    device_states,
    device_parameters,
    static::PSY.StaticInjection,
    dynamic_device::DynamicWrapper{PSY.DynamicGenerator{M, S, A, PSY.DEGOV, P}},
    inner_vars::AbstractVector,
) where {M <: PSY.Machine, S <: PSY.Shaft, A <: PSY.AVR, P <: PSY.PSS}
    tg = PSY.get_prime_mover(dynamic_device)
    #Get mechanical torque to SyncMach
    τm0 = inner_vars[τm_var]
    PSY.set_P_ref!(tg, τm0)
    set_P_ref(dynamic_device, τm0)
    #Update states
    tg_ix = get_local_state_ix(dynamic_device, typeof(tg))
    tg_states = @view device_states[tg_ix]
    tg_states[1] = 0.0
    tg_states[2] = 0.0
    tg_states[3] = 0.0
    tg_states[4] = 0.0
    tg_states[5] = τm0
    return
end

function initialize_tg!(
    device_states,
    p,
    ::PSY.StaticInjection,
    dynamic_device::DynamicWrapper{PSY.DynamicGenerator{M, S, A, PSY.SteamTurbineGov1, P}},
    inner_vars::AbstractVector,
) where {M <: PSY.Machine, S <: PSY.Shaft, A <: PSY.AVR, P <: PSY.PSS}

    #Get mechanical torque to SyncMach
    τm0 = inner_vars[τm_var]
    Δω = 0.0

    tg = PSY.get_prime_mover(dynamic_device)
    #Get Parameters
    params = p[:params][:TurbineGov]
    R = params[:R]
    V_min = params[:valve_position_limits][:min]
    V_max = params[:valve_position_limits][:max]
    inv_R = R < eps() ? 0.0 : (1.0 / R)

    function f!(out, x, params)
        P_ref = x[1]
        x_g1 = x[2]
        x_g2 = x[3]

        R = params[:R]
        T1 = params[:T1]
        V_min = params[:valve_position_limits][:min]
        V_max = params[:valve_position_limits][:max]
        T2 = params[:T2]
        T3 = params[:T3]
        D_T = params[:D_T]
        inv_R = R < eps() ? 0.0 : (1.0 / R)

        ref_in = inv_R * (P_ref - Δω)
        Pm = x_g2 + (T2 / T3) * x_g1

        out[1] = (1.0 / T1) * (ref_in - x_g1) #dx_g1/dt
        x_g1_sat = V_min < x_g1 < V_max ? x_g1 : max(V_min, min(V_max, x_g1))
        out[2] = (1.0 / T3) * (x_g1_sat * (1 - T2 / T3) - x_g2) #dx_g2/dt
        out[3] = (Pm - D_T * Δω) - τm0
    end
    x0 = [1.0 / inv_R, τm0, τm0]
    prob = NonlinearSolve.NonlinearProblem(f!, x0, params)
    sol = NonlinearSolve.solve(
        prob,
        NonlinearSolve.TrustRegion();
        reltol = STRICT_NLSOLVE_F_TOLERANCE,
        abstol = STRICT_NLSOLVE_F_TOLERANCE,
    )
    if !SciMLBase.successful_retcode(sol)
        CRC.@ignore_derivatives @warn("Initialization in TG failed")
    else
        sol_x0 = sol.u
        if (sol_x0[2] >= V_max + BOUNDS_TOLERANCE) ||
           (sol_x0[2] <= V_min - BOUNDS_TOLERANCE)
            CRC.@ignore_derivatives @error(
                "Valve limits for TG in $(PSY.get_name(dynamic_device)) (x_g1 = $(sol_x0[2])), outside its limits V_max = $V_max, Vmin = $V_min.  Consider updating the operating point."
            )
        end
        #Update Control Refs
        PSY.set_P_ref!(tg, sol_x0[1])
        set_P_ref!(p, sol_x0[1])
        #Update states
        tg_ix = get_local_state_ix(dynamic_device, typeof(tg))
        tg_states = @view device_states[tg_ix]
        tg_states[1] = sol_x0[2]
        tg_states[2] = sol_x0[3]
    end
    return
end

function initialize_tg!(
    device_states,
    device_parameters,
    static::PSY.StaticInjection,
    dynamic_device::DynamicWrapper{PSY.DynamicGenerator{M, S, A, PSY.HydroTurbineGov, P}},
    inner_vars::AbstractVector,
) where {M <: PSY.Machine, S <: PSY.Shaft, A <: PSY.AVR, P <: PSY.PSS}

    #Get mechanical torque to SyncMach
    τm0 = inner_vars[τm_var]
    Δω = 0.0

    tg = PSY.get_prime_mover(dynamic_device)
    #Get Parameters
    local_ix_params = get_local_parameter_ix(dynamic_device, PSY.HydroTurbineGov)
    internal_params = @view device_parameters[local_ix_params]
    R,
    r,
    Tr,
    Tf,
    Tg,
    VELM,   #Gate velocity limits not implemented
    G_min,
    G_max,
    Tw,
    At,
    D_T,
    q_nl = internal_params

    function f!(out, x)
        P_ref = x[1]
        x_g1 = x[2]
        x_g2 = x[3]
        x_g3 = x[4]
        x_g4 = x[5]

        c = (1.0 / r) * x_g1 + (1.0 / (r * Tr)) * x_g2
        P_in = P_ref - Δω - R * c
        h = (x_g4 / x_g3)^2

        out[1] = (P_in - x_g1)
        out[2] = x_g1
        out[3] = (c - x_g3)
        out[4] = (1.0 - h)
        out[5] = (x_g4 - q_nl) * h * At - D_T * Δω * x_g3 - τm0
    end
    P0 = R * (q_nl + τm0 / At) # mechanical power initial guess. It migth be different than electrical output power
    x0 = [P0, 0, (r * Tr) * P0 / R, P0 / R, P0 / R]
    sol = NLsolve.nlsolve(f!, x0; ftol = STRICT_NLSOLVE_F_TOLERANCE)
    if !NLsolve.converged(sol)
        CRC.@ignore_derivatives @warn(
            "Initialization of Turbine Governor $(PSY.get_name(static)) failed"
        )
    else
        sol_x0 = sol.zero
        #Error if x_g3 is outside PI limits
        if sol_x0[4] > G_max || sol_x0[4] < G_min
            error(
                "Hydro Turbine Governor $(PSY.get_name(static)) $(sol_x0[4]) outside its gate limits $G_min, $G_max. Consider updating the operating point.",
            )
        end
        #Update Control Refs
        PSY.set_P_ref!(tg, sol_x0[1])
        set_P_ref(dynamic_device, sol_x0[1])
        #Update states
        tg_ix = get_local_state_ix(dynamic_device, typeof(tg))
        tg_states = @view device_states[tg_ix]
        tg_states[1] = sol_x0[2]
        tg_states[2] = sol_x0[3]
        tg_states[3] = sol_x0[4]
        tg_states[4] = sol_x0[5]
    end
    return
end
