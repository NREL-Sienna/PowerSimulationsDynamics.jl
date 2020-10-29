function initialize_tg!(
    device_states,
    static::PSY.StaticInjection,
    dynamic_device::PSY.DynamicGenerator{M, S, A, PSY.TGFixed, P},
) where {M <: PSY.Machine, S <: PSY.Shaft, A <: PSY.AVR, P <: PSY.PSS}
    tg = PSY.get_prime_mover(dynamic_device)
    τm0 = get_inner_vars(dynamic_device)[τm_var]
    eff = PSY.get_efficiency(tg)
    P_ref = τm0 / eff
    PSY.set_P_ref!(tg, P_ref)
    #Update Control Refs
    PSY.get_ext(dynamic_device)[CONTROL_REFS][P_ref_index] = P_ref
end

function initialize_tg!(
    device_states,
    static::PSY.StaticInjection,
    dynamic_device::PSY.DynamicGenerator{M, S, A, PSY.TGTypeI, P},
) where {M <: PSY.Machine, S <: PSY.Shaft, A <: PSY.AVR, P <: PSY.PSS}

    #Get mechanical torque to SyncMach
    τm0 = get_inner_vars(dynamic_device)[τm_var]
    #Get Parameters
    tg = PSY.get_prime_mover(dynamic_device)
    R = PSY.get_R(tg)
    Tc = PSY.get_Tc(tg)
    T3 = PSY.get_T3(tg)
    T4 = PSY.get_T4(tg)
    T5 = PSY.get_T5(tg)

    #Get References
    ω_ref = PSY.get_ext(dynamic_device)[CONTROL_REFS][ω_ref_index]
    ω0 = 1.0

    function f!(out, x)
        P_ref = x[1]
        x_g1 = x[2]
        x_g2 = x[3]
        x_g3 = x[4]

        #Compute auxiliary parameters
        P_in = P_ref + (1.0 / R) * (ω_ref - ω0)

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
        PSY.get_ext(dynamic_device)[CONTROL_REFS][P_ref_index] = sol_x0[1]
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
    static::PSY.StaticInjection,
    dynamic_device::PSY.DynamicGenerator{M, S, A, PSY.TGTypeII, P},
) where {M <: PSY.Machine, S <: PSY.Shaft, A <: PSY.AVR, P <: PSY.PSS}

    #Get mechanical torque to SyncMach
    τm0 = get_inner_vars(dynamic_device)[τm_var]
    #Get parameters
    tg = PSY.get_prime_mover(dynamic_device)
    R = PSY.get_R(tg)
    T1 = PSY.get_T1(tg)
    T2 = PSY.get_T2(tg)
    ω_ref = PSY.get_ext(dynamic_device)[CONTROL_REFS][ω_ref_index]
    ω0 = ω_ref

    function f!(out, x)
        P_ref = x[1]
        xg = x[2]

        out[1] = (1.0 / R) * (T1 / T2) * (ω_ref - ω0) + P_ref / 1.0 + xg - τm0
        out[2] = (1.0 / T2) * ((1.0 / R) * (1 - T2 / T1) * (ω_ref - ω0) - xg)
    end
    x0 = [τm0, 0.0]
    sol = NLsolve.nlsolve(f!, x0)
    if !NLsolve.converged(sol)
        @warn("Initialization in Synch. Machine failed")
    else
        sol_x0 = sol.zero
        #Update Control Refs
        PSY.set_P_ref!(tg, sol_x0[1])
        PSY.get_ext(dynamic_device)[CONTROL_REFS][P_ref_index] = sol_x0[1]
        #Update states
        tg_ix = get_local_state_ix(dynamic_device, PSY.TGTypeII)
        tg_states = @view device_states[tg_ix]
        tg_states[1] = sol_x0[2]
    end
end
