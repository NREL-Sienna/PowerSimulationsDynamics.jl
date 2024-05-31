function mass_matrix_tg_entries!(
    mass_matrix,
    tg::TG,
    global_index::Base.ImmutableDict{Symbol, Int64},
) where {TG <: PSY.TurbineGov}
    @debug "Using default mass matrix entries $TG"
end

function mass_matrix_tg_entries!(
    mass_matrix,
    tg::PSY.HydroTurbineGov,
    global_index::Base.ImmutableDict{Symbol, Int64},
)
    mass_matrix[global_index[:x_g1], global_index[:x_g1]] = PSY.get_Tf(tg)
    mass_matrix[global_index[:x_g3], global_index[:x_g3]] = PSY.get_Tg(tg)
    return
end

function mass_matrix_tg_entries!(
    mass_matrix,
    tg::PSY.DEGOV,
    global_index::Base.ImmutableDict{Symbol, Int64},
)
    mass_matrix[global_index[:x_ecb1], global_index[:x_ecb1]] =
        PSY.get_T1(tg) * PSY.get_T2(tg)
    mass_matrix[global_index[:x_a1], global_index[:x_a1]] = PSY.get_T5(tg) * PSY.get_T6(tg)
    return
end

function mdl_tg_ode!(
    device_states::AbstractArray{<:ACCEPTED_REAL_TYPES},
    ::AbstractArray{<:ACCEPTED_REAL_TYPES},
    p::AbstractArray{<:ACCEPTED_REAL_TYPES},
    inner_vars::AbstractArray{<:ACCEPTED_REAL_TYPES},
    ω_sys::ACCEPTED_REAL_TYPES,
    device::DynamicWrapper{PSY.DynamicGenerator{M, S, A, PSY.TGFixed, P}},
    h,
    t,
) where {M <: PSY.Machine, S <: PSY.Shaft, A <: PSY.AVR, P <: PSY.PSS}
    efficiency = p[:params][:TurbineGov][:efficiency]
    P_ref = p[:refs][:P_ref]
    inner_vars[τm_var] = P_ref * efficiency
    return
end

function mdl_tg_ode!(
    device_states::AbstractArray{<:ACCEPTED_REAL_TYPES},
    output_ode::AbstractArray{<:ACCEPTED_REAL_TYPES},
    device_parameters::AbstractArray{<:ACCEPTED_REAL_TYPES},
    inner_vars::AbstractArray{<:ACCEPTED_REAL_TYPES},
    ω_sys::ACCEPTED_REAL_TYPES,
    device::DynamicWrapper{PSY.DynamicGenerator{M, S, A, PSY.TGTypeI, P}},
    h,
    t,
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
    local_ix_params = get_local_parameter_ix(device, PSY.TGTypeI)
    internal_params = @view device_parameters[local_ix_params]
    R, Ts, Tc, T3, T4, T5, V_min, V_max = internal_params
    inv_R = R < eps() ? 0.0 : (1.0 / R)

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
    device_parameters::AbstractArray{<:ACCEPTED_REAL_TYPES},
    inner_vars::AbstractArray{<:ACCEPTED_REAL_TYPES},
    ω_sys::ACCEPTED_REAL_TYPES,
    device::DynamicWrapper{PSY.DynamicGenerator{M, S, A, PSY.TGTypeII, P}},
    h,
    t,
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

    #Get parameters
    local_ix_params = get_local_parameter_ix(device, PSY.TGTypeII)
    internal_params = @view device_parameters[local_ix_params]
    R, T1, T2 = internal_params
    inv_R = R < eps() ? 0.0 : (1.0 / R)

    #Compute block derivatives
    y_ll1, dxg_dt = lead_lag(inv_R * (ω_ref - ω[1]), xg, 1.0, T1, T2)
    τ_m = y_ll1 + P_ref # Uses P_ref assuming that τ_m0 = P_ref in SS.

    #Compute 1 State TG ODE:
    output_ode[local_ix[1]] = dxg_dt

    #Update mechanical torque
    inner_vars[τm_var] = τ_m

    return
end

function mdl_tg_ode!(
    device_states::AbstractArray{<:ACCEPTED_REAL_TYPES},
    output_ode::AbstractArray{<:ACCEPTED_REAL_TYPES},
    p::AbstractArray{<:ACCEPTED_REAL_TYPES},
    inner_vars::AbstractArray{<:ACCEPTED_REAL_TYPES},
    ω_sys::ACCEPTED_REAL_TYPES,
    device::DynamicWrapper{PSY.DynamicGenerator{M, S, A, PSY.SteamTurbineGov1, P}},
    h,
    t,
) where {M <: PSY.Machine, S <: PSY.Shaft, A <: PSY.AVR, P <: PSY.PSS}

    #Obtain TG
    tg = PSY.get_prime_mover(device)
    #Obtain references
    P_ref = p[:refs][:P_ref]

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
    params = @view(p[:params][:TurbineGov])
    R = params[:R]
    T1 = params[:T1]
    V_min = params[:valve_position_limits][:min]
    V_max = params[:valve_position_limits][:max]
    T2 = params[:T2]
    T3 = params[:T3]
    D_T = params[:D_T]
    inv_R = R < eps() ? 0.0 : (1.0 / R)

    #Compute auxiliary parameters
    ref_in = inv_R * (P_ref - (ω[1] - 1.0))

    #Compute block derivatives
    x_g1_sat, dxg1_dt = low_pass_nonwindup(ref_in, x_g1, 1.0, T1, V_min, V_max)
    y_ll, dxg2_dt = lead_lag(x_g1_sat, x_g2, 1.0, T2, T3)
    P_m = y_ll - D_T * (ω[1] - 1.0)

    #Compute 2 State TG ODE:
    output_ode[local_ix[1]] = dxg1_dt
    output_ode[local_ix[2]] = dxg2_dt

    #Update mechanical torque
    inner_vars[τm_var] = P_m / ω[1]

    return
end

function mdl_tg_ode!(
    device_states::AbstractArray{<:ACCEPTED_REAL_TYPES},
    output_ode::AbstractArray{<:ACCEPTED_REAL_TYPES},
    device_parameters::AbstractArray{<:ACCEPTED_REAL_TYPES},
    inner_vars::AbstractArray{<:ACCEPTED_REAL_TYPES},
    ω_sys::ACCEPTED_REAL_TYPES,
    device::DynamicWrapper{PSY.DynamicGenerator{M, S, A, PSY.GasTG, P}},
    h,
    t,
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
    local_ix_params = get_local_parameter_ix(device, PSY.GasTG)
    internal_params = @view device_parameters[local_ix_params]
    R, T1, T2, T3, AT, Kt, V_min, V_max, D_turb = internal_params
    inv_R = R < eps() ? 0.0 : (1.0 / R)

    #Compute auxiliary parameters
    x_in = min((P_ref - inv_R * (ω[1] - 1.0)), (AT + Kt * (AT - x_g3)))

    #Compute block derivatives
    x_g1_sat, dxg1_dt = low_pass_nonwindup(x_in, x_g1, 1.0, T1, V_min, V_max)
    _, dxg2_dt = low_pass(x_g1_sat, x_g2, 1.0, T2)
    _, dxg3_dt = low_pass(x_g2, x_g3, 1.0, T3)

    #Compute output torque
    P_m = x_g2 - D_turb * (ω[1] - 1.0)

    #Compute 1 State TG ODE:
    output_ode[local_ix[1]] = dxg1_dt
    output_ode[local_ix[2]] = dxg2_dt
    output_ode[local_ix[3]] = dxg3_dt

    #Update mechanical torque
    inner_vars[τm_var] = P_m / ω[1]

    return
end

function mdl_tg_ode!(
    device_states::AbstractArray{<:ACCEPTED_REAL_TYPES},
    output_ode::AbstractArray{<:ACCEPTED_REAL_TYPES},
    device_parameters::AbstractArray{<:ACCEPTED_REAL_TYPES},
    inner_vars::AbstractArray{<:ACCEPTED_REAL_TYPES},
    ω_sys::ACCEPTED_REAL_TYPES,
    device::DynamicWrapper{PSY.DynamicGenerator{M, S, A, PSY.HydroTurbineGov, P}},
    h,
    t,
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
    local_ix_params = get_local_parameter_ix(device, PSY.HydroTurbineGov)
    internal_params = @view device_parameters[local_ix_params]
    R,
    r,
    Tr,
    Tf,
    Tg,
    VELM,
    G_min,
    G_max,
    Tw,
    At,
    D_T,
    q_nl = internal_params

    #Compute block derivatives
    c, dxg2_dt = pi_block_nonwindup(x_g1, x_g2, 1.0 / r, 1.0 / (r * Tr), G_min, G_max)
    P_in = P_ref - Δω - R * c
    _, dxg1_dt = low_pass_mass_matrix(P_in, x_g1, 1.0, Tf)
    _, dxg3_dt = low_pass_mass_matrix(c, x_g3, 1.0, Tg)
    h = (x_g4 / x_g3)^2

    #Compute 4 States TG ODE:
    output_ode[local_ix[1]] = dxg1_dt
    output_ode[local_ix[2]] = clamp(dxg2_dt, -VELM, VELM)
    output_ode[local_ix[3]] = dxg3_dt
    output_ode[local_ix[4]] = (1.0 - h) / Tw

    #Update mechanical torque
    inner_vars[τm_var] = ((x_g4 - q_nl) * h * At - D_T * Δω * x_g3) / ω[1]
    return
end

#Wrapper around h to handle delays of zero
function get_delayed_value(h, t, delay, state, index)
    if delay == 0
        return state
    else
        return h(nothing, t - delay; idxs = index)
    end
end
function mdl_tg_ode!(
    device_states::AbstractArray{<:ACCEPTED_REAL_TYPES},
    output_ode::AbstractArray{<:ACCEPTED_REAL_TYPES},
    device_parameters::AbstractArray{<:ACCEPTED_REAL_TYPES},
    inner_vars::AbstractArray{<:ACCEPTED_REAL_TYPES},
    ω_sys::ACCEPTED_REAL_TYPES,
    device::DynamicWrapper{PSY.DynamicGenerator{M, S, A, PSY.DEGOV, P}},
    h,
    t,
) where {M <: PSY.Machine, S <: PSY.Shaft, A <: PSY.AVR, P <: PSY.PSS}

    #Obtain references
    P_ref = get_P_ref(device)

    #Obtain indices for component w/r to device
    local_ix = get_local_state_ix(device, PSY.DEGOV)

    #Define internal states for component
    internal_states = @view device_states[local_ix]
    x_ecb1 = internal_states[1]
    x_ecb2 = internal_states[2]
    x_a1 = internal_states[3]
    x_a2 = internal_states[4]
    x_a3 = internal_states[5]

    #Obtain external states inputs for component
    external_ix = get_input_port_ix(device, PSY.DEGOV)
    ω = @view device_states[external_ix]

    #Get Parameters
    tg = PSY.get_prime_mover(device)
    T1 = PSY.get_T1(tg)
    T2 = PSY.get_T2(tg)
    T3 = PSY.get_T3(tg)
    K = PSY.get_K(tg)
    T4 = PSY.get_T4(tg)
    T5 = PSY.get_T5(tg)
    T6 = PSY.get_T6(tg)
    Td = PSY.get_Td(tg)

    #Compute block derivatives 
    Δω = ω[1] - 1.0
    y1, dx_ecb1, dx_ecb2 =
        lead_lag_2nd_mass_matrix(-1.0 * Δω, x_ecb1, x_ecb2, T1, T1 * T2, T3, 0.0)
    y2, dx_a1, dx_a2 = lead_lag_2nd_mass_matrix(y1, x_a1, x_a2, T5 + T6, T5 * T6, T4, 0.0)
    y3, dx_a3 = integrator_windup(y2, x_a3, K, 1.0, -Inf, Inf)
    P_m = y3 * (ω[1])
    delayed_x_a3 = get_delayed_value(h, t, Td, x_a3, get_global_index(device)[:x_a3])
    P_m = delayed_x_a3 * (ω[1])

    #Compute 1 State TG ODE:
    output_ode[local_ix[1]] = dx_ecb1
    output_ode[local_ix[2]] = dx_ecb2
    output_ode[local_ix[3]] = dx_a1
    output_ode[local_ix[4]] = dx_a2
    output_ode[local_ix[5]] = dx_a3

    #Update mechanical torque
    inner_vars[τm_var] = P_m / ω[1] #Fails when trying to assign a Dual to a cache of type Float? 

    return
end
