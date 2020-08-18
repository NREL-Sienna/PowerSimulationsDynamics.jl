function mdl_tg_ode!(
    device_states,
    output_ode,
    ω_sys,
    device::PSY.DynamicGenerator{M, S, A, PSY.TGFixed, P},
) where {M <: PSY.Machine, S <: PSY.Shaft, A <: PSY.AVR, P <: PSY.PSS}

    #Update inner vars
    P_ref = PSY.get_ext(device)[CONTROL_REFS][P_ref_index]
    get_inner_vars(device)[τm_var] = P_ref * PSY.get_efficiency(PSY.get_prime_mover(device))

    return
end

function mdl_tg_ode!(
    device_states,
    output_ode,
    ω_sys,
    device::PSY.DynamicGenerator{M, S, A, PSY.TGTypeI, P},
) where {M <: PSY.Machine, S <: PSY.Shaft, A <: PSY.AVR, P <: PSY.PSS}

    #Obtain references
    ω_ref = PSY.get_ext(device)[CONTROL_REFS][ω_ref_index]
    P_ref = PSY.get_ext(device)[CONTROL_REFS][P_ref_index]

    #Obtain indices for component w/r to device
    local_ix = get_local_state_ix(device, PSY.TGTypeI)

    #Define internal states for component
    internal_states = @view device_states[local_ix]
    x_g1 = internal_states[1]
    x_g2 = internal_states[2]
    x_g3 = internal_states[3]

    #Obtain external states inputs for component
    external_ix = get_input_port_ix(device, PSY.TGTypeI)
    ω = @view device_states[external_ix]

    #Get Parameters
    tg = PSY.get_prime_mover(device)
    R = PSY.get_R(tg)
    Ts = PSY.get_Ts(tg)
    Tc = PSY.get_Tc(tg)
    T3 = PSY.get_T3(tg)
    T4 = PSY.get_T4(tg)
    T5 = PSY.get_T5(tg)
    #P_min = PSY.get_P_min(tg)
    #P_max = PSY.get_P_max(tg)

    #Compute auxiliary parameters
    P_in = P_ref + (1.0 / R) * (ω_ref - ω[1])

    #Set anti-windup for P_in. #TODO in callbacks
    #if P_in > P_max
    #    P_in = P_max
    #elseif P_in < P_min
    #    P_in = P_min
    #end

    #Compute 3 States TG ODE:
    output_ode[local_ix[1]] = (1.0 / Ts) * (P_in - x_g1)
    output_ode[local_ix[2]] = (1.0 / Tc) * ((1.0 - T3 / Tc) * x_g1 - x_g2)
    output_ode[local_ix[3]] =
        (1.0 / T5) * ((1.0 - T4 / T5) * (x_g2 + (T3 / Tc) * x_g1) - x_g3)

    #Update mechanical torque
    get_inner_vars(device)[τm_var] = x_g3 + (T4 / T5) * (x_g2 + (T3 / Tc) * x_g1)

    return
end

function mdl_tg_ode!(
    device_states,
    output_ode,
    ω_sys,
    device::PSY.DynamicGenerator{M, S, A, PSY.TGTypeII, P},
) where {M <: PSY.Machine, S <: PSY.Shaft, A <: PSY.AVR, P <: PSY.PSS}

    #Obtain references
    ω_ref = PSY.get_ext(device)[CONTROL_REFS][ω_ref_index]
    P_ref = PSY.get_ext(device)[CONTROL_REFS][P_ref_index]

    #Obtain indices for component w/r to device
    local_ix = get_local_state_ix(device, PSY.TGTypeII)

    #Define internal states for component
    internal_states = @view device_states[local_ix]
    xg = internal_states[1]

    #Obtain external states inputs for component
    external_ix = get_input_port_ix(device, PSY.TGTypeII)
    ω = @view device_states[external_ix]

    #Get Parameters
    tg = PSY.get_prime_mover(device)
    R = PSY.get_R(tg)
    T1 = PSY.get_T1(tg)
    T2 = PSY.get_T2(tg)
    #τ_min = PSY.get_τ_min(tg)
    #τ_max = PSY.get_τ_max(tg)

    #Compute auxiliary parameters
    τ_m = (1.0 / R) * (T1 / T2) * (ω_ref - ω[1]) + P_ref / 1.0 + xg

    #Set anti-windup for τ_m. NOT WORKING
    #if τ_m > τ_max
    #    τ_m = τ_max
    #elseif τ_m < τ_min
    #    τ_m = τ_min
    #end

    #Compute 1 State TG ODE:
    output_ode[local_ix[1]] = (1.0 / T2) * ((1.0 / R) * (1 - T2 / T1) * (ω_ref - ω[1]) - xg)

    #Update mechanical torque
    get_inner_vars(device)[τm_var] = τ_m

    return
end

function mdl_tg_ode!(
    device_states,
    output_ode,
    ω_sys,
    device::PSY.DynamicGenerator{M, S, A, PSY.SteamTurbineGov1, P},
) where {M <: PSY.Machine, S <: PSY.Shaft, A <: PSY.AVR, P <: PSY.PSS}

    #Obtain TG
    tg = PSY.get_prime_mover(device)
    #Obtain references
    P_ref = PSY.get_ext(device)[CONTROL_REFS][P_ref_index]

    #Obtain indices for component w/r to device
    local_ix = get_local_state_ix(device, typeof(tg))

    #Define internal states for component
    internal_states = @view device_states[local_ix]
    x_g1 = internal_states[1]
    x_g2 = internal_states[2]

    #Obtain external states inputs for component
    external_ix = get_input_port_ix(device, typeof(tg))
    ω = @view device_states[external_ix]

    #Get Parameters
    tg = PSY.get_prime_mover(device)
    R = PSY.get_R(tg)
    T1 = PSY.get_T1(tg)
    T2 = PSY.get_T2(tg)
    V_min, V_max = PSY.get_valve_position_limits(tg)
    T3 = PSY.get_T3(tg)
    D_T = PSY.get_D_T(tg)

    #Compute auxiliary parameters
    x_g1_sat = x_g1
    ref_in = (1.0 / R) * (P_ref - (ω[1] - ω_sys))
    Pm = x_g2 + (T2 / T3) * x_g1
    τ_m = Pm - D_T * (ω[1] - ω_sys)

    #Set anti-windup for x_g1
    if x_g1_sat > V_max
        x_g1_sat = V_max
    elseif x_g1_sat < V_min
        x_g1_sat = V_min
    end

    #Compute 2 State TG ODE:
    output_ode[local_ix[1]] = (1.0 / T1) * (ref_in - x_g1) #dx_g1/dt
    output_ode[local_ix[2]] = (1.0 / T3) * (x_g1_sat * (1 - T2 / T3) - x_g2) #dx_g2/dt

    #Update mechanical torque
    get_inner_vars(device)[τm_var] = τ_m

    return
end
