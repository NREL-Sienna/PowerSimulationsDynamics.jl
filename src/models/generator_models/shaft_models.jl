function mass_matrix_shaft_entries!(
    mass_matrix,
    shaft::S,
    global_index::Base.ImmutableDict{Symbol, Int64},
) where {S <: PSY.Shaft}
    @debug "Using default mass matrix entries $S"
end

function mdl_shaft_ode!(
    device_states::AbstractArray{<:ACCEPTED_REAL_TYPES},
    output_ode::AbstractArray{<:ACCEPTED_REAL_TYPES},
    inner_vars::AbstractArray{<:ACCEPTED_REAL_TYPES},
    ω_sys::ACCEPTED_REAL_TYPES,
    dynamic_device::DynamicWrapper{PSY.DynamicGenerator{M, PSY.SingleMass, A, TG, P}},
    h,
    t,
) where {M <: PSY.Machine, A <: PSY.AVR, TG <: PSY.TurbineGov, P <: PSY.PSS}
    f0 = get_system_base_frequency(dynamic_device)
    #Obtain indices for component w/r to device
    local_ix = get_local_state_ix(dynamic_device, PSY.SingleMass)

    #Define internal states for component
    internal_states = @view device_states[local_ix]
    ω = internal_states[2]

    #Obtain inner variables for component
    τe = inner_vars[τe_var]
    τm = inner_vars[τm_var]

    #Get parameters
    shaft = PSY.get_shaft(dynamic_device)
    H = PSY.get_H(shaft)
    D = PSY.get_D(shaft)

    #Compute 2 states ODEs
    output_ode[local_ix[1]] = 2 * π * f0 * (ω - ω_sys)                    #15.5
    output_ode[local_ix[2]] = (1 / (2 * H)) * (τm - τe - D * (ω - 1.0) / ω)   #15.5

    return
end

function mdl_shaft_ode!(
    device_states::AbstractArray{<:ACCEPTED_REAL_TYPES},
    output_ode::AbstractArray{<:ACCEPTED_REAL_TYPES},
    inner_vars::AbstractArray{<:ACCEPTED_REAL_TYPES},
    ω_sys::ACCEPTED_REAL_TYPES,
    dynamic_device::DynamicWrapper{PSY.DynamicGenerator{M, PSY.FiveMassShaft, A, TG, P}},
    h,
    t,
) where {M <: PSY.Machine, A <: PSY.AVR, TG <: PSY.TurbineGov, P <: PSY.PSS}
    f0 = get_system_base_frequency(dynamic_device)

    #Obtain indices for component w/r to device
    local_ix = get_local_state_ix(dynamic_device, PSY.FiveMassShaft)

    #Define internal states for component
    internal_states = @view device_states[local_ix]
    δ = internal_states[1]
    ω = internal_states[2]
    δ_hp = internal_states[3]
    ω_hp = internal_states[4]
    δ_ip = internal_states[5]
    ω_ip = internal_states[6]
    δ_lp = internal_states[7]
    ω_lp = internal_states[8]
    δ_ex = internal_states[9]
    ω_ex = internal_states[10]

    #Obtain inner variables for component
    τe = inner_vars[τe_var]
    τm = inner_vars[τm_var]

    #Get parameters
    shaft = PSY.get_shaft(dynamic_device)
    H = PSY.get_H(shaft)
    H_hp = PSY.get_H_hp(shaft)
    H_ip = PSY.get_H_ip(shaft)
    H_lp = PSY.get_H_lp(shaft)
    H_ex = PSY.get_H_ex(shaft)
    D = PSY.get_D(shaft)
    D_hp = PSY.get_D_hp(shaft)
    D_ip = PSY.get_D_ip(shaft)
    D_lp = PSY.get_D_lp(shaft)
    D_ex = PSY.get_D_ex(shaft)
    D_12 = PSY.get_D_12(shaft)
    D_23 = PSY.get_D_23(shaft)
    D_34 = PSY.get_D_34(shaft)
    D_45 = PSY.get_D_45(shaft)
    K_hp = PSY.get_K_hp(shaft)
    K_ip = PSY.get_K_ip(shaft)
    K_lp = PSY.get_K_lp(shaft)
    K_ex = PSY.get_K_ex(shaft)

    #Compute 10 states ODEs #15.51
    output_ode[local_ix[1]] = 2.0 * π * f0 * (ω - ω_sys)

    output_ode[local_ix[2]] =
        (1.0 / (2.0 * H)) * (
            -τe - D * (ω - 1.0) - D_34 * (ω - ω_lp) - D_45 * (ω - ω_ex) +
            K_lp * (δ_lp - δ) +
            K_ex * (δ_ex - δ)
        )

    output_ode[local_ix[3]] = 2.0 * π * f0 * (ω_hp - ω_sys)

    output_ode[local_ix[4]] =
        (1.0 / (2.0 * H_hp)) *
        (τm - D_hp * (ω_hp - 1.0) - D_12 * (ω_hp - ω_ip) + K_hp * (δ_ip - δ_hp))

    output_ode[local_ix[5]] = 2 * π * f0 * (ω_ip - ω_sys)

    output_ode[local_ix[6]] =
        (1.0 / (2 * H_ip)) * (
            -D_ip * (ω_ip - 1.0) - D_12 * (ω_ip - ω_hp) - D_23 * (ω_ip - ω_lp) +
            K_hp * (δ_hp - δ_ip) +
            K_ip * (δ_lp - δ_ip)
        )

    output_ode[local_ix[7]] = 2.0 * π * f0 * (ω_lp - ω_sys)

    output_ode[local_ix[8]] =
        (1.0 / (2.0 * H_lp)) * (
            -D_lp * (ω_lp - 1.0) - D_23 * (ω_lp - ω_ip) - D_34 * (ω_lp - ω) +
            K_ip * (δ_ip - δ_lp) +
            K_lp * (δ - δ_lp)
        )

    output_ode[local_ix[9]] = 2 * π * f0 * (ω_ex - ω_sys)

    output_ode[local_ix[10]] =
        (1.0 / (2.0 * H_ex)) *
        (-D_ex * (ω_ex - 1.0) - D_45 * (ω_ex - ω) + K_ex * (δ - δ_ex))

    return
end
