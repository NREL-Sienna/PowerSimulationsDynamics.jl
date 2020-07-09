function initialize_tg!(device_states,
    device::PSY.DynamicGenerator{M, S, A, PSY.TGFixed, P},
) where {M <: PSY.Machine, S <: PSY.Shaft, A <: PSY.AVR, P <: PSY.PSS}

    tg = PSY.get_prime_mover(device)
    τm0 = get_inner_vars(device)[τm_var]
    eff = PSY.get_efficiency(tg)
    P_ref = τm0 / eff
    PSY.set_P_ref!(tg, P_ref)
    #Update Control Refs
    device.ext[CONTROL_REFS][P_ref_index] = P_ref
end




function initialize_tg!(device_states,
    device::PSY.DynamicGenerator{M, S, A, PSY.TGTypeII, P},
) where {M <: PSY.Machine, S <: PSY.Shaft, A <: PSY.AVR, P <: PSY.PSS}

    #Get mechanical torque to SyncMach
    τm0 = get_inner_vars(device)[τm_var]
    #Get parameters
    tg = PSY.get_prime_mover(device)
    R = PSY.get_R(tg)
    T1 = PSY.get_T1(tg)
    T2 = PSY.get_T2(tg)
    ω_ref = device.ext[CONTROL_REFS][ω_ref_index]
    ω0 = 1.0

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
        device.ext[CONTROL_REFS][P_ref_index] = sol_x0[1]
        #Update states
        tg_ix = get_local_state_ix(device, PSY.TGTypeII)
        tg_states = @view device_states[tg_ix]
        tg_states[1] = sol_x0[2] 
    end
end