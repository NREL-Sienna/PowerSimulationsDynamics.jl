function mass_matrix_tg_entries!(
    mass_matrix,
    tg::TG,
    global_index::Base.ImmutableDict{Symbol, Int64},
) where {TG <: PSY.TurbineGov}
    @debug "Using default mass matrix entries $TG"
end

function mass_matrix_avr_entries!(
    mass_matrix,
    tg::PSY.HydroTurbineGov,
    global_index::Base.ImmutableDict{Symbol, Int64},
)
    mass_matrix[global_index[:x_g1], global_index[:x_g1]] = PSY.get_Tf(avr)
    mass_matrix[global_index[:x_g3], global_index[:x_g3]] = PSY.get_Tg(avr)
end

function mdl_tg_ode!(
    ::AbstractArray{<:ACCEPTED_REAL_TYPES},
    ::AbstractArray{<:ACCEPTED_REAL_TYPES},
    inner_vars::AbstractArray{<:ACCEPTED_REAL_TYPES},
    ω_sys::ACCEPTED_REAL_TYPES,
    device::DynamicWrapper{PSY.DynamicGenerator{M, S, A, PSY.TGFixed, P}},
) where {M <: PSY.Machine, S <: PSY.Shaft, A <: PSY.AVR, P <: PSY.PSS}

    #Update inner vars
    P_ref = get_P_ref(device)
    inner_vars[τm_var] = P_ref * PSY.get_efficiency(PSY.get_prime_mover(device))

    return
end

function mdl_tg_ode!(
    device_states::AbstractArray{<:ACCEPTED_REAL_TYPES},
    output_ode::AbstractArray{<:ACCEPTED_REAL_TYPES},
    inner_vars::AbstractArray{<:ACCEPTED_REAL_TYPES},
    ω_sys::ACCEPTED_REAL_TYPES,
    device::DynamicWrapper{PSY.DynamicGenerator{M, S, A, PSY.TGTypeI, P}},
) where {M <: PSY.Machine, S <: PSY.Shaft, A <: PSY.AVR, P <: PSY.PSS}

    #Obtain references
    ω_ref = get_ω_ref(device)
    P_ref = get_P_ref(device)

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
    inv_R = PSY.get_R(tg) < eps() ? 0.0 : (1.0 / PSY.get_R(tg))
    Ts = PSY.get_Ts(tg)
    Tc = PSY.get_Tc(tg)
    T3 = PSY.get_T3(tg)
    T4 = PSY.get_T4(tg)
    T5 = PSY.get_T5(tg)
    V_min, V_max = PSY.get_valve_position_limits(tg)

    #Compute auxiliary parameters
    P_in_sat = clamp(P_ref + inv_R * (ω_ref - ω[1]), V_min, V_max)

    #Compute block derivatives
    _, dxg1_dt = low_pass(P_in_sat, x_g1, 1.0, Ts)
    y_ll, dxg2_dt = lead_lag(x_g1, x_g2, 1.0, T3, Tc)
    τ_m, dxg3_dt = lead_lag(y_ll, x_g3, 1.0, T4, T5)

    #Compute 3 States TG ODE:
    output_ode[local_ix[1]] = dxg1_dt
    output_ode[local_ix[2]] = dxg2_dt
    output_ode[local_ix[3]] = dxg3_dt

    #Update mechanical torque
    inner_vars[τm_var] = τ_m

    return
end

function mdl_tg_ode!(
    device_states::AbstractArray{<:ACCEPTED_REAL_TYPES},
    output_ode::AbstractArray{<:ACCEPTED_REAL_TYPES},
    inner_vars::AbstractArray{<:ACCEPTED_REAL_TYPES},
    ω_sys::ACCEPTED_REAL_TYPES,
    device::DynamicWrapper{PSY.DynamicGenerator{M, S, A, PSY.TGTypeII, P}},
) where {M <: PSY.Machine, S <: PSY.Shaft, A <: PSY.AVR, P <: PSY.PSS}

    #Obtain references
    ω_ref = get_ω_ref(device)
    P_ref = get_P_ref(device)

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
    inv_R = PSY.get_R(tg) < eps() ? 0.0 : (1.0 / PSY.get_R(tg))
    T1 = PSY.get_T1(tg)
    T2 = PSY.get_T2(tg)

    #Compute block derivatives
    y_ll1, dxg_dt = lead_lag(inv_R * (ω_ref - ω[1]), xg, 1.0, T1, T2)
    τ_m = y_ll1 + P_ref / 1.0

    #Compute 1 State TG ODE:
    output_ode[local_ix[1]] = dxg_dt

    #Update mechanical torque
    inner_vars[τm_var] = τ_m

    return
end

function mdl_tg_ode!(
    device_states::AbstractArray{<:ACCEPTED_REAL_TYPES},
    output_ode::AbstractArray{<:ACCEPTED_REAL_TYPES},
    inner_vars::AbstractArray{<:ACCEPTED_REAL_TYPES},
    ω_sys::ACCEPTED_REAL_TYPES,
    device::DynamicWrapper{PSY.DynamicGenerator{M, S, A, PSY.SteamTurbineGov1, P}},
) where {M <: PSY.Machine, S <: PSY.Shaft, A <: PSY.AVR, P <: PSY.PSS}

    #Obtain TG
    tg = PSY.get_prime_mover(device)
    #Obtain references
    P_ref = get_P_ref(device)

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
    inv_R = PSY.get_R(tg) < eps() ? 0.0 : (1.0 / PSY.get_R(tg))
    T1 = PSY.get_T1(tg)
    T2 = PSY.get_T2(tg)
    V_min, V_max = PSY.get_valve_position_limits(tg)
    T3 = PSY.get_T3(tg)
    D_T = PSY.get_D_T(tg)

    #Compute auxiliary parameters
    ref_in = inv_R * (P_ref - (ω[1] - 1.0))

    #Compute block derivatives
    x_g1_sat, dxg1_dt = low_pass_nonwindup(ref_in, x_g1, 1.0, T1, V_min, V_max)
    y_ll, dxg2_dt = lead_lag(x_g1_sat, x_g2, 1.0, T2, T3)
    τ_m = y_ll - D_T * (ω[1] - 1.0)

    #Compute 2 State TG ODE:
    output_ode[local_ix[1]] = dxg1_dt
    output_ode[local_ix[2]] = dxg2_dt

    #Update mechanical torque
    inner_vars[τm_var] = τ_m

    return
end

function mdl_tg_ode!(
    device_states::AbstractArray{<:ACCEPTED_REAL_TYPES},
    output_ode::AbstractArray{<:ACCEPTED_REAL_TYPES},
    inner_vars::AbstractArray{<:ACCEPTED_REAL_TYPES},
    ω_sys::ACCEPTED_REAL_TYPES,
    device::DynamicWrapper{PSY.DynamicGenerator{M, S, A, PSY.GasTG, P}},
) where {M <: PSY.Machine, S <: PSY.Shaft, A <: PSY.AVR, P <: PSY.PSS}

    #Obtain references
    P_ref = get_P_ref(device)

    #Obtain indices for component w/r to device
    local_ix = get_local_state_ix(device, PSY.GasTG)

    #Define internal states for component
    internal_states = @view device_states[local_ix]
    x_g1 = internal_states[1]
    x_g2 = internal_states[2]
    x_g3 = internal_states[3]

    #Obtain external states inputs for component
    external_ix = get_input_port_ix(device, PSY.GasTG)
    ω = @view device_states[external_ix]

    #Get Parameters
    tg = PSY.get_prime_mover(device)
    inv_R = PSY.get_R(tg) < eps() ? 0.0 : (1.0 / PSY.get_R(tg))
    T1 = PSY.get_T1(tg)
    T2 = PSY.get_T2(tg)
    T3 = PSY.get_T3(tg)
    D_turb = PSY.get_D_turb(tg)
    AT = PSY.get_AT(tg)
    KT = PSY.get_Kt(tg)
    V_min, V_max = PSY.get_V_lim(tg)

    #Compute auxiliary parameters
    x_in = min((P_ref - inv_R * (ω[1] - 1.0)), (AT + KT * (AT - x_g3)))

    #Compute block derivatives
    x_g1_sat, dxg1_dt = low_pass_nonwindup(x_in, x_g1, 1.0, T1, V_min, V_max)
    _, dxg2_dt = low_pass(x_g1_sat, x_g2, 1.0, T2)
    _, dxg3_dt = low_pass(x_g2, x_g3, 1.0, T3)

    #Compute output torque
    τ_m = x_g2 - D_turb * (ω[1] - 1.0)

    #Compute 1 State TG ODE:
    output_ode[local_ix[1]] = dxg1_dt
    output_ode[local_ix[2]] = dxg2_dt
    output_ode[local_ix[3]] = dxg3_dt

    #Update mechanical torque
    inner_vars[τm_var] = τ_m

    return
end

function mdl_tg_ode!(
    device_states::AbstractArray{<:ACCEPTED_REAL_TYPES},
    output_ode::AbstractArray{<:ACCEPTED_REAL_TYPES},
    inner_vars::AbstractArray{<:ACCEPTED_REAL_TYPES},
    ω_sys::ACCEPTED_REAL_TYPES,
    device::DynamicWrapper{PSY.DynamicGenerator{M, S, A, PSY.HydroTurbineGov, P}},
) where {M <: PSY.Machine, S <: PSY.Shaft, A <: PSY.AVR, P <: PSY.PSS}

    #Obtain references
    P_ref = get_P_ref(device)
    ω_ref = get_ω_ref(device)

    #Obtain indices for component w/r to device
    local_ix = get_local_state_ix(device, PSY.HydroTurbineGov)

    #Define internal states for component
    internal_states = @view device_states[local_ix]
    x_g1 = internal_states[1] # Filter
    x_g2 = internal_states[2] # PI State
    x_g3 = internal_states[3] # Gate
    x_g4 = internal_states[4] # Flow

    #Obtain external states inputs for component
    external_ix = get_input_port_ix(device, PSY.HydroTurbineGov)
    ω = device_states[external_ix][1]
    Δω = ω - ω_ref

    #Get Parameters
    tg = PSY.get_prime_mover(device)
    R = PSY.get_R(tg)
    r = PSY.get_r(tg)
    Tr = PSY.get_Tr(tg)
    Tf = PSY.get_Tf(tg)
    Tg = PSY.get_Tg(tg)
    #Gate velocity limits not implemented
    #VELM = PSY.get_VELM(tg)
    G_min, G_max = PSY.get_gate_position_limits(tg)
    Tw = PSY.get_Tw(tg)
    At = PSY.get_At(tg)
    D_T = PSY.get_D_T(tg)
    q_nl = PSY.get_q_nl(tg)

    #Compute block derivatives
    c, dxg2_dt = pi_block_nonwindup(x_g1, x_g2, 1.0 / r, 1.0 / (r * Tr), G_min, G_max)
    P_in = P_ref - Δω - R * c
    _, dxg1_dt = low_pass_mass_matrix(P_in, x_g1, 1.0, Tf)
    _, dxg3_dt = low_pass_mass_matrix(c, x_g3, 1.0, Tg)
    h = (x_g4 / x_g3)^2

    #Compute 4 States TG ODE:
    output_ode[local_ix[1]] = dxg1_dt
    output_ode[local_ix[2]] = dxg2_dt
    output_ode[local_ix[3]] = dxg3_dt
    output_ode[local_ix[4]] = (1.0 - h) / Tw

    #Update mechanical torque
    inner_vars[τm_var] = (x_g4 - q_nl) * h * At - D_T * Δω * x_g3
    return
end
