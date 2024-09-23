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
    p,
    static::PSY.StaticInjection,
    dynamic_device::DynamicWrapper{PSY.DynamicGenerator{M, S, A, PSY.TGTypeI, P}},
    inner_vars::AbstractVector,
) where {M <: PSY.Machine, S <: PSY.Shaft, A <: PSY.AVR, P <: PSY.PSS}

    #Get mechanical torque to SyncMach
    τm0 = inner_vars[τm_var]

    tg = PSY.get_prime_mover(dynamic_device)
    #Get Parameters
    params = p[:params][:TurbineGov]
    R = params[:R]
    Tc = params[:Tc]
    T3 = params[:T3]
    T4 = params[:T4]
    T5 = params[:T5]
    V_min = params[:valve_position_limits][:min]
    V_max = params[:valve_position_limits][:max]

    #Get References
    ω_ref = p[:refs][:ω_ref]
    ω0 = 1.0
    function f!(out, x, params)
        P_ref = x[1]
        x_g1 = x[2]
        x_g2 = x[3]
        x_g3 = x[4]

        R = params[:R]
        Tc = params[:Tc]
        T3 = params[:T3]
        T4 = params[:T4]
        T5 = params[:T5]
        ω0 = 1.0
        #Compute auxiliary parameters
        inv_R = R < eps() ? 0.0 : (1.0 / R)
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
    prob = NonlinearSolve.NonlinearProblem{true}(f!, x0, params)
    sol = NonlinearSolve.solve(
        prob,
        NonlinearSolve.TrustRegion();
        sensealg = SciMLSensitivity.SteadyStateAdjoint(),
        reltol = STRICT_NLSOLVE_F_TOLERANCE,
        abstol = STRICT_NLSOLVE_F_TOLERANCE,
    )
    if !SciMLBase.successful_retcode(sol)
        @warn("Initialization of Turbine Governor $(PSY.get_name(static)) failed")
    else
        sol_x0 = sol.u
        #Update Control Refs
        PSY.set_P_ref!(tg, sol_x0[1])
        set_P_ref!(p, sol_x0[1])
        #Update states
        tg_ix = get_local_state_ix(dynamic_device, PSY.TGTypeI)
        tg_states = @view device_states[tg_ix]
        if (sol_x0[2] > V_max) || (sol_x0[2] < V_min)
            @error(
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
    p,
    static::PSY.StaticInjection,
    dynamic_device::DynamicWrapper{PSY.DynamicGenerator{M, S, A, PSY.TGTypeII, P}},
    inner_vars::AbstractVector,
) where {M <: PSY.Machine, S <: PSY.Shaft, A <: PSY.AVR, P <: PSY.PSS}

    #Get mechanical torque to SyncMach
    τm0 = inner_vars[τm_var]
    tg = PSY.get_prime_mover(dynamic_device)

    params = [
        p[:params][:TurbineGov][:R],
        p[:params][:TurbineGov][:T1],
        p[:params][:TurbineGov][:T2],
        p[:refs][:ω_ref],
    ]

    function f!(out, x, params)
        P_ref = x[1]
        xg = x[2]
        R, T1, T2, ω_ref = params
        ω0 = ω_ref
        inv_R = R < eps() ? 0.0 : (1.0 / R)
        out[1] = inv_R * (T1 / T2) * (ω_ref - ω0) + P_ref / 1.0 + xg - τm0
        out[2] = (1.0 / T2) * (inv_R * (1 - T2 / T1) * (ω_ref - ω0) - xg)
    end
    x0 = [τm0, 0.0]
    prob = NonlinearSolve.NonlinearProblem{true}(f!, x0, params)
    sol = NonlinearSolve.solve(
        prob,
        NonlinearSolve.TrustRegion();
        sensealg = SciMLSensitivity.SteadyStateAdjoint(),
        reltol = STRICT_NLSOLVE_F_TOLERANCE,
        abstol = STRICT_NLSOLVE_F_TOLERANCE,
    )
    if !SciMLBase.successful_retcode(sol)
        @warn("Initialization of Turbine Governor $(PSY.get_name(static)) failed")
    else
        sol_x0 = sol.u
        #Update Control Refs
        PSY.set_P_ref!(tg, sol_x0[1])
        set_P_ref!(p, sol_x0[1])
        #Update states
        tg_ix = get_local_state_ix(dynamic_device, PSY.TGTypeII)
        tg_states = @view device_states[tg_ix]
        tg_states[1] = sol_x0[2]
    end
    return
end

function initialize_tg!(
    device_states,
    p,
    static::PSY.StaticInjection,
    dynamic_device::DynamicWrapper{PSY.DynamicGenerator{M, S, A, PSY.GasTG, P}},
    inner_vars::AbstractVector,
) where {M <: PSY.Machine, S <: PSY.Shaft, A <: PSY.AVR, P <: PSY.PSS}

    #Get mechanical torque to SyncMach
    τm0 = inner_vars[τm_var]
    Δω = 0.0

    tg = PSY.get_prime_mover(dynamic_device)
    #Get parameters
    params = p[:params][:TurbineGov]
    R = params[:R]
    inv_R = R < eps() ? 0.0 : (1.0 / R)

    function f!(out, x, params)
        P_ref = x[1]
        x_g1 = x[2]
        x_g2 = x[3]
        x_g3 = x[4]

        AT = params[:AT]
        Kt = params[:Kt]
        V_min = params[:V_lim][:min]
        V_max = params[:V_lim][:max]
        D_turb = params[:D_turb]

        x_in = min((P_ref - inv_R * Δω), (AT + Kt * (AT - x_g3)))
        out[1] = (x_in - x_g1) #dx_g1/dt
        x_g1_sat = V_min < x_g1 < V_max ? x_g1 : max(V_min, min(V_max, x_g1))
        out[2] = (x_g1_sat - x_g2)
        out[3] = (x_g2 - x_g3)
        out[4] = (x_g2 - D_turb * Δω) - τm0
    end
    x0 = [1.0 / inv_R, τm0, τm0, τm0]
    prob = NonlinearSolve.NonlinearProblem{true}(f!, x0, params)
    sol = NonlinearSolve.solve(
        prob,
        NonlinearSolve.TrustRegion();
        sensealg = SciMLSensitivity.SteadyStateAdjoint(),
        reltol = STRICT_NLSOLVE_F_TOLERANCE,
        abstol = STRICT_NLSOLVE_F_TOLERANCE,
    )
    if !SciMLBase.successful_retcode(sol)
        @warn("Initialization in TG failed")
    else
        sol_x0 = sol.u
        #Update Control Refs
        PSY.set_P_ref!(tg, sol_x0[1])
        set_P_ref!(p, sol_x0[1])
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
    p,
    static::PSY.StaticInjection,
    dynamic_device::DynamicWrapper{PSY.DynamicGenerator{M, S, A, PSY.DEGOV, P}},
    inner_vars::AbstractVector,
) where {M <: PSY.Machine, S <: PSY.Shaft, A <: PSY.AVR, P <: PSY.PSS}
    tg = PSY.get_prime_mover(dynamic_device)
    #Get mechanical torque to SyncMach
    τm0 = inner_vars[τm_var]
    PSY.set_P_ref!(tg, τm0)
    set_P_ref!(p, τm0)
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
    prob = NonlinearSolve.NonlinearProblem{true}(f!, x0, params)
    sol = NonlinearSolve.solve(
        prob,
        NonlinearSolve.TrustRegion();
        sensealg = SciMLSensitivity.SteadyStateAdjoint(),
        reltol = STRICT_NLSOLVE_F_TOLERANCE,
        abstol = STRICT_NLSOLVE_F_TOLERANCE,
    )
    if !SciMLBase.successful_retcode(sol)
        @warn("Initialization in TG failed")
    else
        sol_x0 = sol.u
        if (sol_x0[2] >= V_max + BOUNDS_TOLERANCE) ||
           (sol_x0[2] <= V_min - BOUNDS_TOLERANCE)
            @error(
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
    p,
    static::PSY.StaticInjection,
    dynamic_device::DynamicWrapper{PSY.DynamicGenerator{M, S, A, PSY.HydroTurbineGov, P}},
    inner_vars::AbstractVector,
) where {M <: PSY.Machine, S <: PSY.Shaft, A <: PSY.AVR, P <: PSY.PSS}

    #Get mechanical torque to SyncMach
    τm0 = inner_vars[τm_var]
    Δω = 0.0

    tg = PSY.get_prime_mover(dynamic_device)
    #Get Parameters
    params = p[:params][:TurbineGov]
    R = params[:R]
    r = params[:r]
    Tr = params[:Tr]
    G_min = params[:gate_position_limits][:min]
    G_max = params[:gate_position_limits][:max]
    At = params[:At]
    q_nl = params[:q_nl]

    function f!(out, x, params)
        P_ref = x[1]
        x_g1 = x[2]
        x_g2 = x[3]
        x_g3 = x[4]
        x_g4 = x[5]

        R = params[:R]
        r = params[:r]
        Tr = params[:Tr]
        At = params[:At]
        D_T = params[:D_T]
        q_nl = params[:q_nl]
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
    prob = NonlinearSolve.NonlinearProblem{true}(f!, x0, params)
    sol = NonlinearSolve.solve(
        prob,
        NonlinearSolve.TrustRegion();
        sensealg = SciMLSensitivity.SteadyStateAdjoint(),
        reltol = STRICT_NLSOLVE_F_TOLERANCE,
        abstol = STRICT_NLSOLVE_F_TOLERANCE,
    )
    if !SciMLBase.successful_retcode(sol)
        @warn("Initialization of Turbine Governor $(PSY.get_name(static)) failed")
    else
        sol_x0 = sol.u
        #Error if x_g3 is outside PI limits
        if sol_x0[4] > G_max || sol_x0[4] < G_min
            error(
                "Hydro Turbine Governor $(PSY.get_name(static)) $(sol_x0[4]) outside its gate limits $G_min, $G_max. Consider updating the operating point.",
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
        tg_states[3] = sol_x0[4]
        tg_states[4] = sol_x0[5]
    end
    return
end

function initialize_tg!(
    device_states,
    static::PSY.StaticInjection,
    dynamic_device::DynamicWrapper{PSY.DynamicGenerator{M, S, A, PSY.PIDGOV, P}},
    inner_vars::AbstractVector,
) where {M <: PSY.Machine, S <: PSY.Shaft, A <: PSY.AVR, P <: PSY.PSS}

    #Get mechanical torque to SyncMach
    τm0 = inner_vars[τm_var]
    P0 = PSY.get_active_power(static)
    Δω = 0.0
    #Get parameters
    tg = PSY.get_prime_mover(dynamic_device)
    feedback_flag = PSY.get_feedback_flag(tg)
    Rperm = PSY.get_Rperm(tg)
    Kp = PSY.get_Kp(tg)
    Ki = PSY.get_Ki(tg)
    Kd = PSY.get_Kd(tg)
    Ta = PSY.get_Ta(tg)
    Tb = PSY.get_Tb(tg)
    gate_openings = PSY.get_gate_openings(tg)
    power_gate_openings = PSY.get_power_gate_openings(tg)
    G_min, G_max = PSY.get_G_lim(tg)
    A_tw = PSY.get_A_tw(tg)
    Tw = PSY.get_Tw(tg)

    function f!(out, x)
        P_ref = x[1]
        x_g1 = x[2]
        x_g2 = x[3]
        x_g3 = x[4]
        x_g4 = x[5]
        x_g5 = x[6]
        x_g6 = x[7]
        x_g7 = x[8]

        if feedback_flag == 0
            P_in = P_ref - P0
        else
            P_in = P_ref - x_g6
        end
        pid_input = x_g1
        pi_out, dxg2_dt = pi_block(pid_input, x_g2, Kp, Ki)
        pd_out, dxg4_dt = high_pass_mass_matrix(pid_input, x_g4, Kd, Ta)

        integrator_input = (1.0 / Tb) * (x_g5 - x_g6)
        power_at_gate =
            three_level_gate_to_power_map(x_g6, gate_openings, power_gate_openings)
        Tz = A_tw * Tw
        y_LL_out, dxg7_dt = lead_lag(power_at_gate, x_g7, 1.0, -Tz, Tz / 2.0)

        out[1] = y_LL_out - τm0
        out[2] = (Rperm * P_in - x_g1)
        out[3] = dxg2_dt
        out[4] = pi_out - x_g3
        out[5] = dxg4_dt
        out[6] = x_g3 + pd_out - x_g5
        out[7] = integrator_input
        out[8] = dxg7_dt
    end
    gate0 = three_level_power_to_gate_map(τm0, gate_openings, power_gate_openings)
    if feedback_flag == 0
        P_ref_guess = P0
    else
        P_ref_guess = gate0
    end
    x0 = [P_ref_guess, 0.0, gate0, gate0, 0.0, gate0, gate0, 3.0 * τm0]
    sol = NLsolve.nlsolve(f!, x0; ftol = STRICT_NLSOLVE_F_TOLERANCE)
    if !NLsolve.converged(sol)
        @warn("Initialization of Turbine Governor $(PSY.get_name(static)) failed")
    else
        sol_x0 = sol.zero
        #Error if x_g3 is outside PI limits
        if sol_x0[7] > G_max || sol_x0[7] < G_min
            @error(
                "PIDGOV Turbine Governor $(PSY.get_name(static)) $(sol_x0[4]) outside its gate limits $G_min, $G_max. Consider updating the operating point.",
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
        tg_states[5] = sol_x0[6]
        tg_states[6] = sol_x0[7]
        tg_states[7] = sol_x0[8]
    end
    return
end

function initialize_tg!(
    device_states,
    static::PSY.StaticInjection,
    dynamic_device::DynamicWrapper{PSY.DynamicGenerator{M, S, A, PSY.WPIDHY, P}},
    inner_vars::AbstractVector,
) where {M <: PSY.Machine, S <: PSY.Shaft, A <: PSY.AVR, P <: PSY.PSS}

    #Get mechanical torque to SyncMach
    τm0 = inner_vars[τm_var]
    τe0 = inner_vars[τe_var]
    P0 = PSY.get_active_power(static)

    #Get parameters
    tg = PSY.get_prime_mover(dynamic_device)
    reg = PSY.get_reg(tg)
    Kp = PSY.get_Kp(tg)
    Ki = PSY.get_Ki(tg)
    Kd = PSY.get_Kd(tg)
    Ta = PSY.get_Ta(tg)
    gate_openings = PSY.get_gate_openings(tg)
    power_gate_openings = PSY.get_power_gate_openings(tg)
    G_min, G_max = PSY.get_G_lim(tg)
    P_min, P_max = PSY.get_P_lim(tg)
    Tw = PSY.get_Tw(tg)

    #Compute controller parameters for equivalent TF
    Kp_prime = (-Ta * Ki) + Kp
    Kd_prime = (Ta^2 * Ki) - (Ta * Kp) + Kd
    Ki_prime = Ki

    function f!(out, x)
        P_ref = x[1]
        x_g1 = x[2]
        x_g2 = x[3]
        x_g3 = x[4]
        x_g4 = x[5]
        x_g5 = x[6]
        x_g6 = x[7]
        x_g7 = x[8]

        pid_input = x_g1
        pi_out, dxg2_dt = pi_block(pid_input, x_g2, Kp_prime, Ki_prime)
        pd_out, dxg4_dt = high_pass(pid_input, x_g4, Kd_prime, Ta)

        power_at_gate =
            three_level_gate_to_power_map(x_g6, gate_openings, power_gate_openings)

        y_LL_out, dxg7_dt = lead_lag(power_at_gate, x_g7, 1.0, -Tw, Tw / 2.0)

        out[1] = y_LL_out - τm0
        out[2] = (P0 - P_ref) * reg - x_g1
        out[3] = dxg2_dt
        out[4] = pi_out + pd_out - x_g3
        out[5] = dxg4_dt
        out[6] = x_g3 - x_g5
        out[7] = x_g5
        out[8] = dxg7_dt
    end
    gate0 = three_level_power_to_gate_map(τm0, gate_openings, power_gate_openings)
    P_ref_guess = P0
    xg1_guess = τe0 - P_ref_guess

    #x0 = [P_ref_guess, xg1_guess, 0.0, 0.0, 0.0, 0.0, gate0, 3.0 * τm0]
    x0 = [P_ref_guess, 0.0, gate0, gate0, 0.0, gate0, gate0, 3.0 * τm0]
    sol = NLsolve.nlsolve(f!, x0; ftol = STRICT_NLSOLVE_F_TOLERANCE)
    if !NLsolve.converged(sol)
        @warn("Initialization of Turbine Governor $(PSY.get_name(static)) failed")
    else
        sol_x0 = sol.zero
        #Error if x_g3 is outside PI limits
        if sol_x0[7] > G_max || sol_x0[7] < G_min
            @error(
                "WPIDHY Turbine Governor $(PSY.get_name(static)) $(sol_x0[7]) outside its gate limits $G_min, $G_max. Consider updating the operating point.",
            )
        elseif τm0 > P_max || τm0 < P_min
            @error(
                "WPIDHY Turbine Governor $(PSY.get_name(static)) $(sol_x0[8]) outside its power limits $P_min, $P_max. Consider updating the operating point.",
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
        tg_states[5] = sol_x0[6]
        tg_states[6] = sol_x0[7]
        tg_states[7] = sol_x0[8]
    end
    return
end
