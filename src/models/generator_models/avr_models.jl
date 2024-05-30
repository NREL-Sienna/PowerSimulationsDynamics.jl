function mass_matrix_avr_entries!(
    mass_matrix,
    avr::T,
    global_index::Base.ImmutableDict{Symbol, Int64},
) where {T <: PSY.AVR}
    CRC.@ignore_derivatives @debug "Using default mass matrix entries $T"
    return
end

function mass_matrix_avr_entries!(
    mass_matrix,
    avr::PSY.SEXS,
    global_index::Base.ImmutableDict{Symbol, Int64},
)
    mass_matrix[global_index[:Vf], global_index[:Vf]] = PSY.get_Te(avr)
    mass_matrix[global_index[:Vr], global_index[:Vr]] = PSY.get_Tb(avr)
    return
end

function mass_matrix_avr_entries!(
    mass_matrix,
    avr::PSY.SCRX,
    global_index::Base.ImmutableDict{Symbol, Int64},
)
    mass_matrix[global_index[:Vr1], global_index[:Vr1]] = PSY.get_Tb(avr) # left hand side
    mass_matrix[global_index[:Vr2], global_index[:Vr2]] = PSY.get_Te(avr) #
    return
end

function mass_matrix_avr_entries!(
    mass_matrix,
    avr::PSY.EXST1,
    global_index::Base.ImmutableDict{Symbol, Int64},
)
    mass_matrix[global_index[:Vm], global_index[:Vm]] = PSY.get_Tr(avr)
    mass_matrix[global_index[:Vrll], global_index[:Vrll]] = PSY.get_Tb(avr)
    mass_matrix[global_index[:Vr], global_index[:Vr]] = PSY.get_Ta(avr)
    return
end

function mass_matrix_avr_entries!(
    mass_matrix,
    avr::PSY.EXAC1,
    global_index::Base.ImmutableDict{Symbol, Int64},
)
    mass_matrix[global_index[:Vm], global_index[:Vm]] = PSY.get_Tr(avr)
    mass_matrix[global_index[:Vr1], global_index[:Vr1]] = PSY.get_Tb(avr)
    mass_matrix[global_index[:Vr2], global_index[:Vr2]] = PSY.get_Ta(avr)
    return
end

function mass_matrix_avr_entries!(
    mass_matrix,
    avr::PSY.ESST1A,
    global_index::Base.ImmutableDict{Symbol, Int64},
)
    mass_matrix[global_index[:Vm], global_index[:Vm]] = PSY.get_Tr(avr)
    mass_matrix[global_index[:Vr1], global_index[:Vr1]] = PSY.get_Tb(avr)
    mass_matrix[global_index[:Vr2], global_index[:Vr2]] = PSY.get_Tb1(avr)
    mass_matrix[global_index[:Va], global_index[:Va]] = PSY.get_Ta(avr)
    return
end

function mdl_avr_ode!(
    ::AbstractArray{<:ACCEPTED_REAL_TYPES},
    ::AbstractArray{<:ACCEPTED_REAL_TYPES},
    p::AbstractArray{<:ACCEPTED_REAL_TYPES},
    inner_vars::AbstractArray{<:ACCEPTED_REAL_TYPES},
    dynamic_device::DynamicWrapper{PSY.DynamicGenerator{M, S, PSY.AVRFixed, TG, P}},
    h,
    t,
) where {M <: PSY.Machine, S <: PSY.Shaft, TG <: PSY.TurbineGov, P <: PSY.PSS}
    V_ref = p[:refs][:V_ref]
    #Update Vf voltage on inner vars. In AVRFixed, Vf = V_ref
    inner_vars[Vf_var] = V_ref
    return
end

function mdl_avr_ode!(
    device_states::AbstractArray{<:ACCEPTED_REAL_TYPES},
    output_ode::AbstractArray{<:ACCEPTED_REAL_TYPES},
    device_parameters::AbstractArray{<:ACCEPTED_REAL_TYPES},
    inner_vars::AbstractArray{<:ACCEPTED_REAL_TYPES},
    dynamic_device::DynamicWrapper{PSY.DynamicGenerator{M, S, PSY.AVRSimple, TG, P}},
    h,
    t,
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
    local_ix_params = get_local_parameter_ix(dynamic_device, PSY.AVRSimple)
    internal_params = @view device_parameters[local_ix_params]
    Kv = internal_params[1]

    #Compute ODEs
    output_ode[local_ix[1]] = Kv * (V_ref - V_th)

    #Update inner_vars
    inner_vars[Vf_var] = Vf

    return
end

function mdl_avr_ode!(
    device_states::AbstractArray{<:ACCEPTED_REAL_TYPES},
    output_ode::AbstractArray{<:ACCEPTED_REAL_TYPES},
    device_parameters::AbstractArray{<:ACCEPTED_REAL_TYPES},
    inner_vars::AbstractArray{<:ACCEPTED_REAL_TYPES},
    dynamic_device::DynamicWrapper{PSY.DynamicGenerator{M, S, PSY.AVRTypeI, TG, P}},
    h,
    t,
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
    local_ix_params = get_local_parameter_ix(dynamic_device, PSY.AVRTypeI)
    internal_params = @view device_parameters[local_ix_params]
    Ka, Ke, Kf, Ta, Te, Tf, Tr, Ae, Be = internal_params

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
    device_parameters::AbstractArray{<:ACCEPTED_REAL_TYPES},
    inner_vars,
    dynamic_device::DynamicWrapper{PSY.DynamicGenerator{M, S, PSY.AVRTypeII, TG, P}},
    h,
    t,
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
    local_ix_params = get_local_parameter_ix(dynamic_device, PSY.AVRTypeII)
    internal_params = @view device_parameters[local_ix_params]
    K0,
    T1,
    T2,
    T3,
    T4,
    Te,
    Tr,
    Va_min,
    Va_max,
    Ae,
    Be = internal_params

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
    device_parameters::AbstractArray{<:ACCEPTED_REAL_TYPES},
    inner_vars::AbstractArray{<:ACCEPTED_REAL_TYPES},
    dynamic_device::DynamicWrapper{PSY.DynamicGenerator{M, S, PSY.ESAC1A, TG, P}},
    h,
    t,
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
    local_ix_params = get_local_parameter_ix(dynamic_device, PSY.ESAC1A)
    internal_params = @view device_parameters[local_ix_params]
    Tr,
    Tb,
    Tc,
    Ka,
    Ta,
    Va_min,
    Va_max,
    Te,
    Kf,
    Tf,
    Kc,
    Kd,
    Ke,
    Vr_min,
    Vr_max = internal_params
    inv_Tr = Tr < eps() ? 1.0 : 1.0 / Tr
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
    p::AbstractArray{<:ACCEPTED_REAL_TYPES},
    inner_vars::AbstractArray{<:ACCEPTED_REAL_TYPES},
    dynamic_device::DynamicWrapper{PSY.DynamicGenerator{M, S, PSY.SEXS, TG, P}},
    h,
    t,
) where {M <: PSY.Machine, S <: PSY.Shaft, TG <: PSY.TurbineGov, P <: PSY.PSS}
    #Obtain references
    V0_ref = p[:refs][:V_ref]
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
    params = @view(p[:params][:AVR])
    Ta_Tb = params[:Ta_Tb]
    Tb = params[:Tb]
    K = params[:K]
    Te = params[:Te]
    V_min = params[:V_lim][:min]
    V_max = params[:V_lim][:max]
    Ta = Tb * Ta_Tb

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

function mdl_avr_ode!(
    device_states::AbstractArray{<:ACCEPTED_REAL_TYPES},
    output_ode::AbstractArray{<:ACCEPTED_REAL_TYPES},
    device_parameters::AbstractArray{<:ACCEPTED_REAL_TYPES},
    inner_vars::AbstractArray{<:ACCEPTED_REAL_TYPES},
    dynamic_device::DynamicWrapper{PSY.DynamicGenerator{M, S, PSY.SCRX, TG, P}},
    h,
    t,
) where {M <: PSY.Machine, S <: PSY.Shaft, TG <: PSY.TurbineGov, P <: PSY.PSS}

    #Obtain references
    V0_ref = get_V_ref(dynamic_device)

    #Obtain indices for component w/r to device
    local_ix = get_local_state_ix(dynamic_device, PSY.SCRX) #

    #Define inner states for component
    internal_states = @view device_states[local_ix]
    Vr1 = internal_states[1]
    Vr2 = internal_states[2]

    #Define external states for device
    V_th = sqrt(inner_vars[VR_gen_var]^2 + inner_vars[VI_gen_var]^2) # real and imaginary >> Ec
    Vs = inner_vars[V_pss_var] #Vs, PSS output
    Ifd = inner_vars[Xad_Ifd_var] # read Lad Ifd

    #Get parameters << keep
    avr = PSY.get_avr(dynamic_device)
    Ta_Tb = PSY.get_Ta_Tb(avr) # Ta/Tb
    Tb = PSY.get_Tb(avr)
    Ta = Tb * Ta_Tb
    Te = PSY.get_Te(avr)
    K = PSY.get_K(avr)
    V_min, V_max = PSY.get_Efd_lim(avr) #Efd_lim (V_lim)
    switch = PSY.get_switch(avr) # reads switch parameters
    rc_rfd = PSY.get_rc_rfd(avr)

    #Compute auxiliary parameters << keep
    V_in = V0_ref + Vs - V_th #sum of V
    V_LL, dVr1_dt = lead_lag_mass_matrix(V_in, Vr1, 1.0, Ta, Tb) # 1st block
    Vr2_sat, dVr2_dt = low_pass_nonwindup_mass_matrix(V_LL, Vr2, K, Te, V_min, V_max) # gain K , 2nd block

    # Switch multiplier
    mult = switch == 0 ? V_th : one(typeof(V_th))
    V_ex = mult * Vr2_sat

    #Negative current logic
    if rc_rfd == 0.0 # a float
        E_fd = V_ex
    else
        E_fd = Ifd > 0.0 ? V_ex : -Ifd * rc_rfd
    end

    #Compute 2 States AVR ODE: << move this after? (final computation)
    output_ode[local_ix[1]] = dVr1_dt
    output_ode[local_ix[2]] = dVr2_dt

    #Update inner_vars << do this after
    inner_vars[Vf_var] = E_fd # field voltage from rc_rfd

    return
end

function mdl_avr_ode!(
    device_states::AbstractArray,
    output_ode::AbstractArray,
    device_parameters::AbstractArray{<:ACCEPTED_REAL_TYPES},
    inner_vars::AbstractArray,
    dynamic_device::DynamicWrapper{PSY.DynamicGenerator{M, S, PSY.EXST1, TG, P}},
    h,
    t,
) where {M <: PSY.Machine, S <: PSY.Shaft, TG <: PSY.TurbineGov, P <: PSY.PSS}

    #Obtain references
    V0_ref = get_V_ref(dynamic_device)

    #Obtain indices for component w/r to device
    local_ix = get_local_state_ix(dynamic_device, PSY.EXST1)

    #Define inner states for component
    internal_states = @view device_states[local_ix]
    Vm = internal_states[1]
    Vrll = internal_states[2]
    Vr = internal_states[3]
    Vfb = internal_states[4]

    #Define external states for device
    Vt = sqrt(inner_vars[VR_gen_var]^2 + inner_vars[VI_gen_var]^2) # machine's terminal voltage
    Vs = inner_vars[V_pss_var] # PSS output 
    Ifd = inner_vars[Xad_Ifd_var] # machine's field current in exciter base 

    #Get Parameters
    local_ix_params = get_local_parameter_ix(dynamic_device, PSY.EXST1)
    internal_params = @view device_parameters[local_ix_params]
    Tr,
    Vi_min,
    Vi_max,
    Tc,
    Tb,
    Ka,
    Ta,
    Vr_min,
    Vr_max,
    Kc,
    Kf,
    Tf = internal_params

    #Compute auxiliary parameters
    V_ref = V0_ref + Vs

    # Compute block derivatives
    _, dVm_dt = low_pass_mass_matrix(Vt, Vm, 1.0, Tr)
    y_hp, dVfb_dt = high_pass(Vr, Vfb, Kf, Tf)
    y_ll, dVrll_dt =
        lead_lag_mass_matrix(clamp(V_ref - Vm - y_hp, Vi_min, Vi_max), Vrll, 1.0, Tc, Tb)
    _, dVr_dt = low_pass_mass_matrix(y_ll, Vr, Ka, Ta)

    #Compute 4 States AVR ODE:
    output_ode[local_ix[1]] = dVm_dt
    output_ode[local_ix[2]] = dVrll_dt
    output_ode[local_ix[3]] = dVr_dt
    output_ode[local_ix[4]] = dVfb_dt

    #Update inner_vars
    Vf = clamp(Vr, Vt * Vr_min - Kc * Ifd, Vt * Vr_max - Kc * Ifd)
    inner_vars[Vf_var] = Vf
    return
end

function mdl_avr_ode!(
    device_states::AbstractArray,
    output_ode::AbstractArray,
    device_parameters::AbstractArray{<:ACCEPTED_REAL_TYPES},
    inner_vars::AbstractArray,
    dynamic_device::DynamicWrapper{PSY.DynamicGenerator{M, S, PSY.EXAC1, TG, P}},
    h,
    t,
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
    V_th = sqrt(inner_vars[VR_gen_var]^2 + inner_vars[VI_gen_var]^2) # machine's terminal voltage
    Vs = inner_vars[V_pss_var] # PSS output 
    Xad_Ifd = inner_vars[Xad_Ifd_var] # machine's field current in exciter base

    #Get parameters
    local_ix_params = get_local_parameter_ix(dynamic_device, PSY.EXAC1)
    internal_params = @view device_parameters[local_ix_params]
    Tr,
    Tb,
    Tc,
    Ka,
    Ta,
    Vr_min,
    Vr_max,
    Te,
    Kf,
    Tf,
    Kc,
    Kd,
    Ke = internal_params

    #Obtain saturation
    Se = saturation_function(avr, Ve)

    #Compute auxiliary parameters
    I_N = Kc * Xad_Ifd / Ve
    V_FE = Kd * Xad_Ifd + Ke * Ve + Se * Ve
    Vf = Ve * rectifier_function(I_N)

    #Compute block derivatives
    _, dVm_dt = low_pass_mass_matrix(V_th, Vm, 1.0, Tr)
    V_F, dVr3_dt = high_pass(V_FE, Vr3, Kf, Tf)
    V_in = V_ref + Vs - Vm - V_F
    y_ll, dVr1_dt = lead_lag_mass_matrix(V_in, Vr1, 1.0, Tc, Tb)
    y_Vr, dVr2_dt = low_pass_nonwindup_mass_matrix(y_ll, Vr2, Ka, Ta, Vr_min, Vr_max)
    _, dVe_dt = integrator_nonwindup(y_Vr - V_FE, Ve, 1.0, Te, 0.0, Inf)

    #Compute 4 States AVR ODE:
    output_ode[local_ix[1]] = dVm_dt
    output_ode[local_ix[2]] = dVr1_dt
    output_ode[local_ix[3]] = dVr2_dt
    output_ode[local_ix[4]] = dVe_dt
    output_ode[local_ix[5]] = dVr3_dt

    #Update inner_vars
    inner_vars[Vf_var] = Vf
    return
end

function mdl_avr_ode!(
    device_states::AbstractArray,
    output_ode::AbstractArray,
    device_parameters::AbstractArray{<:ACCEPTED_REAL_TYPES},
    inner_vars::AbstractArray,
    dynamic_device::DynamicWrapper{PSY.DynamicGenerator{M, S, PSY.ESST1A, TG, P}},
    h,
    t,
) where {M <: PSY.Machine, S <: PSY.Shaft, TG <: PSY.TurbineGov, P <: PSY.PSS}

    #Obtain references
    V0_ref = get_V_ref(dynamic_device)

    #Obtain indices for component w/r to device
    local_ix = get_local_state_ix(dynamic_device, PSY.ESST1A)

    #Define inner states for component
    internal_states = @view device_states[local_ix]
    Vm = internal_states[1]
    Vr1 = internal_states[2]
    Vr2 = internal_states[3]
    Va = internal_states[4]
    Vr3 = internal_states[5]

    #Define external states for device
    Vt = sqrt(inner_vars[VR_gen_var]^2 + inner_vars[VI_gen_var]^2) # machine's terminal voltage
    Vs = inner_vars[V_pss_var] # PSS output 
    Ifd = inner_vars[Xad_Ifd_var] # machine's field current in exciter base 

    #Get parameters
    avr = PSY.get_avr(dynamic_device)
    UEL = PSY.get_UEL_flags(avr)
    VOS = PSY.get_PSS_flags(avr)
    Tr = PSY.get_Tr(avr)
    Vi_min, Vi_max = PSY.get_Vi_lim(avr)
    Tc = PSY.get_Tc(avr)
    Tb = PSY.get_Tb(avr)
    Tc1 = PSY.get_Tc1(avr)
    Tb1 = PSY.get_Tb1(avr)
    Ka = PSY.get_Ka(avr)
    Ta = PSY.get_Ta(avr)
    Va_min, Va_max = PSY.get_Va_lim(avr)
    Vr_min, Vr_max = PSY.get_Vr_lim(avr)
    Kc = PSY.get_Kc(avr)
    Kf = PSY.get_Kf(avr)
    Tf = PSY.get_Tf(avr)
    K_lr = PSY.get_K_lr(avr)
    I_lr = PSY.get_I_lr(avr)

    #Compute auxiliary parameters
    Itemp = K_lr * (Ifd - I_lr)
    Iresult = Itemp > 0.0 ? Itemp : 0.0

    if VOS == 1
        V_ref = V0_ref + Vs
        Va_sum = Va - Iresult
    elseif VOS == 2
        V_ref = V0_ref
        Va_sum = Va - Iresult + Vs
    end

    # Compute block derivatives
    _, dVm_dt = low_pass_mass_matrix(Vt, Vm, 1.0, Tr)
    y_hp, dVr3_dt = high_pass(Va_sum, Vr3, Kf, Tf)
    y_ll1, dVr1_dt =
        lead_lag_mass_matrix(clamp(V_ref - Vm - y_hp, Vi_min, Vi_max), Vr1, 1.0, Tc, Tb)
    y_ll2, dVr2_dt =
        lead_lag_mass_matrix(y_ll1, Vr2, 1.0, Tc1, Tb1)
    _, dVa_dt = low_pass_nonwindup_mass_matrix(y_ll2, Va, Ka, Ta, Va_min, Va_max)

    #Compute 5 States AVR ODE:
    output_ode[local_ix[1]] = dVm_dt
    output_ode[local_ix[2]] = dVr1_dt
    output_ode[local_ix[3]] = dVr2_dt
    output_ode[local_ix[4]] = dVa_dt
    output_ode[local_ix[5]] = dVr3_dt

    #Update inner_vars
    Vf = clamp(Va_sum, Vt * Vr_min, Vt * Vr_max - Kc * Ifd)
    inner_vars[Vf_var] = Vf
    return
end
