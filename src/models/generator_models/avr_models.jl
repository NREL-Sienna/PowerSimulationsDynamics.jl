function mass_matrix_avr_entries!(
    mass_matrix,
    avr::T,
    global_index::Base.ImmutableDict{Symbol, Int64},
) where {T <: PSY.AVR}
    @debug "Using default mass matrix entries $T"
end

function mass_matrix_avr_entries!(
    mass_matrix,
    avr::PSY.SEXS,
    global_index::Base.ImmutableDict{Symbol, Int64},
)
    mass_matrix[global_index[:Vf], global_index[:Vf]] = PSY.get_Te(avr)
    mass_matrix[global_index[:Vr], global_index[:Vr]] = PSY.get_Tb(avr)
end

function mdl_avr_ode!(
    ::AbstractArray{<:ACCEPTED_REAL_TYPES},
    ::AbstractArray{<:ACCEPTED_REAL_TYPES},
    inner_vars::AbstractArray{<:ACCEPTED_REAL_TYPES},
    dynamic_device::DynamicWrapper{PSY.DynamicGenerator{M, S, PSY.AVRFixed, TG, P}},
) where {M <: PSY.Machine, S <: PSY.Shaft, TG <: PSY.TurbineGov, P <: PSY.PSS}

    #Update Vf voltage on inner vars. In AVRFixed, Vf = V_ref
    inner_vars[Vf_var] = get_V_ref(dynamic_device)

    return
end

function mdl_avr_ode!(
    device_states::AbstractArray{<:ACCEPTED_REAL_TYPES},
    output_ode::AbstractArray{<:ACCEPTED_REAL_TYPES},
    inner_vars::AbstractArray{<:ACCEPTED_REAL_TYPES},
    dynamic_device::DynamicWrapper{PSY.DynamicGenerator{M, S, PSY.AVRSimple, TG, P}},
) where {M <: PSY.Machine, S <: PSY.Shaft, TG <: PSY.TurbineGov, P <: PSY.PSS}

    #Obtain references
    V_ref = get_V_ref(dynamic_device)

    #Obtain indices for component w/r to device
    local_ix = get_local_state_ix(dynamic_device, PSY.AVRSimple)

    #Define inner states for component
    internal_states = @view device_states[local_ix]
    Vf = internal_states[1]

    #Define external states for device
    V_th = sqrt(inner_vars[VR_gen_var]^2 + inner_vars[VI_gen_var]^2)

    #Get Parameters
    Kv = PSY.get_Kv(PSY.get_avr(dynamic_device))

    #Compute ODEs
    output_ode[local_ix[1]] = Kv * (V_ref - V_th)

    #Update inner_vars
    inner_vars[Vf_var] = Vf

    return
end

function mdl_avr_ode!(
    device_states::AbstractArray{<:ACCEPTED_REAL_TYPES},
    output_ode::AbstractArray{<:ACCEPTED_REAL_TYPES},
    inner_vars::AbstractArray{<:ACCEPTED_REAL_TYPES},
    dynamic_device::DynamicWrapper{PSY.DynamicGenerator{M, S, PSY.AVRTypeI, TG, P}},
) where {M <: PSY.Machine, S <: PSY.Shaft, TG <: PSY.TurbineGov, P <: PSY.PSS}

    #Obtain references
    V0_ref = get_V_ref(dynamic_device)

    #Obtain indices for component w/r to device
    local_ix = get_local_state_ix(dynamic_device, PSY.AVRTypeI)

    #Define inner states for component
    internal_states = @view device_states[local_ix]
    Vf = internal_states[1]
    Vr1 = internal_states[2]
    Vr2 = internal_states[3]
    Vm = internal_states[4]

    #Define external states for device
    V_th = sqrt(inner_vars[VR_gen_var]^2 + inner_vars[VI_gen_var]^2)
    Vs = inner_vars[V_pss_var]

    #Get parameters
    avr = PSY.get_avr(dynamic_device)
    Ka = PSY.get_Ka(avr)
    Ke = PSY.get_Ke(avr)
    Kf = PSY.get_Kf(avr)
    Ta = PSY.get_Ta(avr)
    Te = PSY.get_Te(avr)
    Tf = PSY.get_Tf(avr)
    Tr = PSY.get_Tr(avr)
    Ae = PSY.get_Ae(avr)
    Be = PSY.get_Be(avr)

    #Compute auxiliary parameters
    Se_Vf = Ae * exp(Be * abs(Vf)) #16.13
    V_ref = V0_ref + Vs

    # Compute block derivatives
    _, dVm_dt = low_pass(V_th, Vm, 1.0, Tr)
    y_hp, dVr2_dt = high_pass(Vf, Vr2, Kf, Tf)
    _, dVr1_dt = low_pass(V_ref - Vm - y_hp, Vr1, Ka, Ta)
    _, dVf_dt = low_pass_modified(Vr1, Vf, 1.0, Ke + Se_Vf, Te)

    #Compute 4 States AVR ODE:
    output_ode[local_ix[1]] = dVf_dt #16.12c
    output_ode[local_ix[2]] = dVr1_dt #16.12a
    output_ode[local_ix[3]] = dVr2_dt #16.12b
    output_ode[local_ix[4]] = dVm_dt #16.11

    #Update inner_vars
    inner_vars[Vf_var] = Vf

    return
end

function mdl_avr_ode!(
    device_states::AbstractArray{<:ACCEPTED_REAL_TYPES},
    output_ode::AbstractArray{<:ACCEPTED_REAL_TYPES},
    inner_vars,
    dynamic_device::DynamicWrapper{PSY.DynamicGenerator{M, S, PSY.AVRTypeII, TG, P}},
) where {M <: PSY.Machine, S <: PSY.Shaft, TG <: PSY.TurbineGov, P <: PSY.PSS}

    #Obtain references
    V0_ref = get_V_ref(dynamic_device)

    #Obtain indices for component w/r to device
    local_ix = get_local_state_ix(dynamic_device, PSY.AVRTypeII)

    #Define inner states for component
    internal_states = @view device_states[local_ix]
    Vf = internal_states[1]
    Vr1 = internal_states[2]
    Vr2 = internal_states[3]
    Vm = internal_states[4]

    #Define external states for device
    V_th = sqrt(inner_vars[VR_gen_var]^2 + inner_vars[VI_gen_var]^2)
    Vs = inner_vars[V_pss_var]

    #Get parameters
    avr = PSY.get_avr(dynamic_device)
    K0 = PSY.get_K0(avr)
    T1 = PSY.get_T1(avr)
    T2 = PSY.get_T2(avr)
    T3 = PSY.get_T3(avr)
    T4 = PSY.get_T4(avr)
    Te = PSY.get_Te(avr)
    Tr = PSY.get_Tr(avr)
    Ae = PSY.get_Ae(avr)
    Be = PSY.get_Be(avr)
    Va_min, Va_max = PSY.get_Va_lim(avr)

    #Compute auxiliary parameters
    Se_Vf = Ae * exp(Be * abs(Vf)) #16.13
    V_ref = V0_ref + Vs

    # Compute block derivatives
    _, dVm_dt = low_pass(V_th, Vm, 1.0, Tr)
    y_ll1, dVr1_dt = lead_lag(V_ref - Vm, Vr1, K0, T2, T1)
    y_ll2, dVr2_dt = lead_lag(y_ll1, K0 * Vr2, 1.0, K0 * T4, K0 * T3)
    Vr = clamp(y_ll2, Va_min, Va_max)
    _, dVf_dt = low_pass_modified(Vr, Vf, 1.0, 1.0 + Se_Vf, Te)

    #Compute 4 States AVR ODE:
    output_ode[local_ix[1]] = dVf_dt #16.18
    output_ode[local_ix[2]] = dVr1_dt #16.14
    output_ode[local_ix[3]] = dVr2_dt
    output_ode[local_ix[4]] = dVm_dt #16.11

    #Update inner_vars
    inner_vars[Vf_var] = Vf

    return
end

function mdl_avr_ode!(
    device_states::AbstractArray{<:ACCEPTED_REAL_TYPES},
    output_ode::AbstractArray{<:ACCEPTED_REAL_TYPES},
    inner_vars::AbstractArray{<:ACCEPTED_REAL_TYPES},
    dynamic_device::DynamicWrapper{PSY.DynamicGenerator{M, S, PSY.ESAC1A, TG, P}},
) where {M <: PSY.Machine, S <: PSY.Shaft, TG <: PSY.TurbineGov, P <: PSY.PSS}

    #Obtain references
    V_ref = get_V_ref(dynamic_device)

    #Obtain avr
    avr = PSY.get_avr(dynamic_device)

    #Obtain indices for component w/r to device
    local_ix = get_local_state_ix(dynamic_device, typeof(avr))

    #Define inner states for component
    internal_states = @view device_states[local_ix]
    Vm = internal_states[1]
    Vr1 = internal_states[2]
    Vr2 = internal_states[3]
    Ve = internal_states[4]
    Vr3 = internal_states[5]

    #Define external states for device
    V_th = sqrt(inner_vars[VR_gen_var]^2 + inner_vars[VI_gen_var]^2)
    Vs = inner_vars[V_pss_var]
    Xad_Ifd = inner_vars[Xad_Ifd_var]

    #Get parameters
    Tr = PSY.get_Tr(avr)
    Ta = PSY.get_Ta(avr)
    Tb = PSY.get_Tb(avr)
    Tc = PSY.get_Tc(avr)
    inv_Tr = Tr < eps() ? 1.0 : 1.0 / Tr
    Ka = PSY.get_Ka(avr)
    Va_min, Va_max = PSY.get_Va_lim(avr) #Not used without UEL or OEL
    Te = PSY.get_Te(avr) # Te > 0
    Kf = PSY.get_Kf(avr)
    Tf = PSY.get_Tf(avr) # Te > 0
    Kc = PSY.get_Kc(avr)
    Kd = PSY.get_Kd(avr)
    Ke = PSY.get_Ke(avr)
    Vr_min, Vr_max = PSY.get_Vr_lim(avr)
    #Obtain saturation
    Se = saturation_function(avr, Ve)

    #Compute auxiliary parameters
    I_N = Kc * Xad_Ifd / Ve
    V_FE = Kd * Xad_Ifd + Ke * Ve + Se * Ve
    Vf = Ve * rectifier_function(I_N)

    # Compute blocks
    _, dVm_dt = low_pass(V_th, Vm, 1.0, 1.0 / inv_Tr)
    V_F, dVr3_dt = high_pass(V_FE, Vr3, Kf, Tf)
    V_in = V_ref + Vs - Vm - V_F
    y_ll, dVr1_dt = lead_lag(V_in, Vr1, 1.0, Tc, Tb)
    _, dVr2_dt = low_pass_nonwindup(y_ll, Vr2, Ka, Ta, Va_min, Va_max)

    #Set clamping for Vr2.
    V_R = clamp(Vr2, Vr_min, Vr_max)

    #Compute 4 States AVR ODE:
    output_ode[local_ix[1]] = dVm_dt
    output_ode[local_ix[2]] = dVr1_dt #dVr1/dt
    output_ode[local_ix[3]] = dVr2_dt
    output_ode[local_ix[4]] = (1.0 / Te) * (V_R - V_FE) #dVe/dt
    output_ode[local_ix[5]] = dVr3_dt

    #Update inner_vars
    inner_vars[Vf_var] = Vf

    return
end

function mdl_avr_ode!(
    device_states::AbstractArray{<:ACCEPTED_REAL_TYPES},
    output_ode::AbstractArray{<:ACCEPTED_REAL_TYPES},
    inner_vars::AbstractArray{<:ACCEPTED_REAL_TYPES},
    dynamic_device::DynamicWrapper{PSY.DynamicGenerator{M, S, PSY.SEXS, TG, P}},
) where {M <: PSY.Machine, S <: PSY.Shaft, TG <: PSY.TurbineGov, P <: PSY.PSS}

    #Obtain references
    V0_ref = get_V_ref(dynamic_device)

    #Obtain indices for component w/r to device
    local_ix = get_local_state_ix(dynamic_device, PSY.SEXS)

    #Define inner states for component
    internal_states = @view device_states[local_ix]
    Vf = internal_states[1]
    Vr = internal_states[2]

    #Define external states for device
    V_th = sqrt(inner_vars[VR_gen_var]^2 + inner_vars[VI_gen_var]^2)
    Vs = inner_vars[V_pss_var]

    #Get parameters
    avr = PSY.get_avr(dynamic_device)
    Ta_Tb = PSY.get_Ta_Tb(avr)
    Tb = PSY.get_Tb(avr)
    Ta = Tb * Ta_Tb
    Te = PSY.get_Te(avr)
    K = PSY.get_K(avr)
    V_min, V_max = PSY.get_V_lim(avr)

    #Compute auxiliary parameters
    V_in = V0_ref + Vs - V_th
    V_LL, dVr_dt = lead_lag_mass_matrix(V_in, Vr, 1.0, Ta, Tb)
    _, dVf_dt = low_pass_nonwindup_mass_matrix(V_LL, Vf, K, Te, V_min, V_max)

    #Compute 2 States AVR ODE:
    output_ode[local_ix[1]] = dVf_dt
    output_ode[local_ix[2]] = dVr_dt

    #Update inner_vars
    inner_vars[Vf_var] = Vf

    return
end
