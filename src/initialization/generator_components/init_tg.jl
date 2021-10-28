function initialize_tg!(
    device_states,
    ::PSY.StaticInjection,
    dynamic_device::DynamicWrapper{PSY.DynamicGenerator{M, S, A, PSY.TGFixed, P}},
    inner_vars::AbstractVector,
) where {M <: PSY.Machine, S <: PSY.Shaft, A <: PSY.AVR, P <: PSY.PSS}
    tg = PSY.get_prime_mover(dynamic_device)
    τm0 = inner_vars[τm_var]
    eff = PSY.get_efficiency(tg)
    P_ref = τm0 / eff
    #PSY.set_P_ref!(tg, P_ref)
    #Update Control Refs
    set_P_ref(dynamic_device, P_ref)
end

function initialize_tg!(
    device_states,
    ::PSY.StaticInjection,
    dynamic_device::DynamicWrapper{PSY.DynamicGenerator{M, S, A, PSY.TGTypeI, P}},
    inner_vars::AbstractVector,
) where {M <: PSY.Machine, S <: PSY.Shaft, A <: PSY.AVR, P <: PSY.PSS}

    #Get mechanical torque to SyncMach
    τm0 = inner_vars[τm_var]
    #Get Parameters
    tg = PSY.get_prime_mover(dynamic_device)
    inv_R = PSY.get_R(tg) < eps() ? 0.0 : (1.0 / PSY.get_R(tg))
    Tc = PSY.get_Tc(tg)
    T3 = PSY.get_T3(tg)
    T4 = PSY.get_T4(tg)
    T5 = PSY.get_T5(tg)

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
        out[3] = ((1.0 - T4 / T5) * (x_g2 + (T3 / Tc) * x_g1) - x_g3)
        out[4] = x_g3 + (T4 / T5) * (x_g2 + (T3 / Tc) * x_g1) - τm0
    end
    x0 = [τm0, τm0, (1.0 - T3 / Tc) * τm0, τm0]
    sol = NLsolve.nlsolve(f!, x0)
    if !NLsolve.converged(sol)
        @warn("Initialization in Synch. Machine failed")
    else
        sol_x0 = sol.zero
        #Update Control Refs
        PSY.set_P_ref!(tg, sol_x0[1])
        set_P_ref(dynamic_device, sol_x0[1])
        #Update states
        tg_ix = get_local_state_ix(dynamic_device, PSY.TGTypeI)
        tg_states = @view device_states[tg_ix]
        tg_states[1] = sol_x0[2]
        tg_states[2] = sol_x0[3]
        tg_states[3] = sol_x0[4]
    end
end

function initialize_tg!(
    device_states,
    ::PSY.StaticInjection,
    dynamic_device::DynamicWrapper{PSY.DynamicGenerator{M, S, A, PSY.TGTypeII, P}},
    inner_vars::AbstractVector,
) where {M <: PSY.Machine, S <: PSY.Shaft, A <: PSY.AVR, P <: PSY.PSS}

    #Get mechanical torque to SyncMach
    τm0 = inner_vars[τm_var]
    #Get parameters
    tg = PSY.get_prime_mover(dynamic_device)
    inv_R = PSY.get_R(tg) < eps() ? 0.0 : (1.0 / PSY.get_R(tg))
    T1 = PSY.get_T1(tg)
    T2 = PSY.get_T2(tg)
    ω_ref = get_ω_ref(dynamic_device)
    ω0 = ω_ref

    function f!(out, x)
        P_ref = x[1]
        xg = x[2]
        out[1] = inv_R * (T1 / T2) * (ω_ref - ω0) + P_ref / 1.0 + xg - τm0
        out[2] = (1.0 / T2) * (inv_R * (1 - T2 / T1) * (ω_ref - ω0) - xg)
    end
    x0 = [τm0, 0.0]
    sol = NLsolve.nlsolve(f!, x0)
    if !NLsolve.converged(sol)
        @warn("Initialization in Synch. Machine failed")
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
end

function initialize_tg!(
    device_states,
    ::PSY.StaticInjection,
    dynamic_device::DynamicWrapper{PSY.DynamicGenerator{M, S, A, PSY.GasTG, P}},
    inner_vars::AbstractVector,
) where {M <: PSY.Machine, S <: PSY.Shaft, A <: PSY.AVR, P <: PSY.PSS}

    #Get mechanical torque to SyncMach
    τm0 = inner_vars[τm_var]
    Δω = 0.0
    #Get parameters
    tg = PSY.get_prime_mover(dynamic_device)
    inv_R = PSY.get_R(tg) < eps() ? 0.0 : (1.0 / PSY.get_R(tg))
    T1 = PSY.get_T1(tg)
    T2 = PSY.get_T2(tg)
    T3 = PSY.get_T3(tg)
    D_turb = PSY.get_D_turb(tg)
    AT = PSY.get_AT(tg)
    KT = PSY.get_Kt(tg)
    V_min, V_max = PSY.get_V_lim(tg)
    ω_ref = get_ω_ref(dynamic_device)
    ω0 = ω_ref

    function f!(out, x)
        P_ref = x[1]
        x_g1 = x[2]
        x_g2 = x[3]
        x_g3 = x[4]

        x_in = min((P_ref - inv_R * Δω), (AT + KT * (AT - x_g3)))
        out[1] = (x_in - x_g1) #dx_g1/dt
        x_g1_sat = V_min < x_g1 < V_max ? x_g1 : max(V_min, min(V_max, x_g1))
        out[2] = (x_g1_sat - x_g2)
        out[3] = (x_g2 - x_g3)
        out[4] = (x_g2 - D_turb * Δω) - τm0
    end
    x0 = [1.0 / inv_R, τm0, τm0, τm0]
    sol = NLsolve.nlsolve(f!, x0)
    if !NLsolve.converged(sol)
        @warn("Initialization in Synch. Machine failed")
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
end

function initialize_tg!(
    device_states,
    ::PSY.StaticInjection,
    dynamic_device::DynamicWrapper{PSY.DynamicGenerator{M, S, A, PSY.SteamTurbineGov1, P}},
    inner_vars::AbstractVector,
) where {M <: PSY.Machine, S <: PSY.Shaft, A <: PSY.AVR, P <: PSY.PSS}

    #Get mechanical torque to SyncMach
    τm0 = inner_vars[τm_var]
    Δω = 0.0
    #Get parameters
    tg = PSY.get_prime_mover(dynamic_device)
    inv_R = PSY.get_R(tg) < eps() ? 0.0 : (1.0 / PSY.get_R(tg))
    T1 = PSY.get_T1(tg)
    T2 = PSY.get_T2(tg)
    T3 = PSY.get_T3(tg)
    V_min, V_max = PSY.get_valve_position_limits(tg)
    D_T = PSY.get_D_T(tg)

    function f!(out, x)
        P_ref = x[1]
        x_g1 = x[2]
        x_g2 = x[3]

        ref_in = inv_R * (P_ref - Δω)
        Pm = x_g2 + (T2 / T3) * x_g1

        out[1] = (1.0 / T1) * (ref_in - x_g1) #dx_g1/dt
        x_g1_sat = V_min < x_g1 < V_max ? x_g1 : max(V_min, min(V_max, x_g1))
        out[2] = (1.0 / T3) * (x_g1_sat * (1 - T2 / T3) - x_g2) #dx_g2/dt
        out[3] = (Pm - D_T * Δω) - τm0
    end
    x0 = [1.0 / inv_R, τm0, τm0]
    sol = NLsolve.nlsolve(f!, x0)
    if !NLsolve.converged(sol)
        @warn("Initialization in TG failed")
    else
        sol_x0 = sol.zero
        if (sol_x0[2] > V_max) || (sol_x0[2] < V_min)
            @error(
                "Valve limits for TG in $(PSY.get_name(dynamic_device)) are bounded (x_g1 = $(sol_x0[2])). Consider updating their values."
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
    end
end

function initialize_tg!(
    device_states,
    static::PSY.StaticInjection,
    dynamic_device::DynamicWrapper{PSY.DynamicGenerator{M, S, A, PSY.HydroTurbineGov, P}},
    inner_vars::AbstractVector,
) where {M <: PSY.Machine, S <: PSY.Shaft, A <: PSY.AVR, P <: PSY.PSS}

    #Get mechanical torque to SyncMach
    τm0 = inner_vars[τm_var]
    Δω = 0.0
    #Get parameters
    tg = PSY.get_prime_mover(dynamic_device)
    R = PSY.get_R(tg)
    r = PSY.get_r(tg)
    Tr = PSY.get_Tr(tg)
    #Gate velocity limits not implemented
    #VELM = PSY.get_VELM(tg)
    G_min, G_max = PSY.get_gate_position_limits(tg)
    Tw = PSY.get_Tw(tg)
    At = PSY.get_At(tg)
    D_T = PSY.get_D_T(tg)
    q_nl = PSY.get_q_nl(tg)

    function f!(out, x)
        P_ref = x[1]
        x_g1 = x[2]
        x_g2 = x[3]
        x_g3 = x[4]
        x_g4 = x[5]

        c = (1.0 / r) * x_g1 + (1.0 / (r * Tr)) * x_g2
        P_in = P_ref - Δω - R * c
        h = (x_g4 / x_g3)^2

        out[1] = (P_in - x_g1) #dx_g1/dt
        out[2] = x_g1
        out[3] = c - x_g3
        out[4] = (1.0 - h) / Tw
        out[5] = (x_g4 - q_nl) * h * At - D_T * Δω * x_g3 - τm0
    end
    P0 = PSY.get_active_power(static)
    x0 = [P0, 0.0, P0 / R, P0 / (R * r * Tr), τm0]
    sol = NLsolve.nlsolve(f!, x0)
    if !NLsolve.converged(sol)
        @warn("Initialization in Synch. Machine failed")
    else
        sol_x0 = sol.zero
        #Error if x_g3 is outside PI limits
        if sol_x0[4] > G_max || sol_x0[4] < G_min
            error(
                "Hydro Turbine Governor of device $(PSY.get_name(static)) is outside its gate limits. Check parameters or Power Flow",
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
end
