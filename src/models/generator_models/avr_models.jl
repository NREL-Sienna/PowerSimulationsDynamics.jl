##################################
###### Mass Matrix Entries #######
##################################

function mass_matrix_avr_entries!(
    mass_matrix,
    avr::T,
    global_index::Base.ImmutableDict{Symbol, Int64},
) where {T <: PSY.AVR}
    @debug "Using default mass matrix entries $T"
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

function mass_matrix_avr_entries!(
    mass_matrix,
    avr::PSY.ST6B,
    global_index::Base.ImmutableDict{Symbol, Int64},
)
    mass_matrix[global_index[:Vm], global_index[:Vm]] = PSY.get_Tr(avr)
    mass_matrix[global_index[:x_d], global_index[:x_d]] = PSY.get_T_da(avr)
    return
end

function mass_matrix_avr_entries!(
    mass_matrix,
    avr::PSY.ST8C,
    global_index::Base.ImmutableDict{Symbol, Int64},
)
    mass_matrix[global_index[:Vm], global_index[:Vm]] = PSY.get_Tr(avr)
    mass_matrix[global_index[:x_a3], global_index[:x_a3]] = PSY.get_T_a(avr)
    mass_matrix[global_index[:x_a4], global_index[:x_a4]] = PSY.get_T_f(avr)
    return
end

##################################
##### Differential Equations #####
##################################

function mdl_avr_ode!(
    ::AbstractArray{<:ACCEPTED_REAL_TYPES},
    ::AbstractArray{<:ACCEPTED_REAL_TYPES},
    inner_vars::AbstractArray{<:ACCEPTED_REAL_TYPES},
    dynamic_device::DynamicWrapper{PSY.DynamicGenerator{M, S, PSY.AVRFixed, TG, P}},
    h,
    t,
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
    h,
    t,
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

function mdl_avr_ode!(
    device_states::AbstractArray{<:ACCEPTED_REAL_TYPES},
    output_ode::AbstractArray{<:ACCEPTED_REAL_TYPES},
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

    #Get parameters
    avr = PSY.get_avr(dynamic_device)
    Tr = PSY.get_Tr(avr)
    Vi_min, Vi_max = PSY.get_Vi_lim(avr)
    Tc = PSY.get_Tc(avr)
    Tb = PSY.get_Tb(avr)
    Ka = PSY.get_Ka(avr)
    Ta = PSY.get_Ta(avr)
    Vr_min, Vr_max = PSY.get_Vr_lim(avr)
    Kc = PSY.get_Kc(avr)
    Kf = PSY.get_Kf(avr)
    Tf = PSY.get_Tf(avr)

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
    Tr = PSY.get_Tr(avr)
    Tb = PSY.get_Tb(avr)
    Tc = PSY.get_Tc(avr)
    Ka = PSY.get_Ka(avr)
    Ta = PSY.get_Ta(avr)
    Vr_min, Vr_max = PSY.get_Vr_lim(avr)
    Te = PSY.get_Te(avr) # Te > 0
    Kf = PSY.get_Kf(avr)
    Tf = PSY.get_Tf(avr) # Tf > 0
    Kc = PSY.get_Kc(avr)
    Kd = PSY.get_Kd(avr)
    Ke = PSY.get_Ke(avr)

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

######################################################################
function mdl_avr_ode!(
    device_states::AbstractArray,
    output_ode::AbstractArray,
    inner_vars::AbstractArray,
    dynamic_device::DynamicWrapper{PSY.DynamicGenerator{M, S, PSY.ST6B, TG, P}},
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
    x_i = internal_states[2]
    x_d = internal_states[3]
    Vg = internal_states[4]

    #Define external states for device
    V_th = sqrt(inner_vars[VR_gen_var]^2 + inner_vars[VI_gen_var]^2) # machine's terminal voltage
    Vs = inner_vars[V_pss_var] # PSS output 
    Xad_Ifd = inner_vars[Xad_Ifd_var] # machine's field current in exciter base

    #Get parameters
    Tr = PSY.get_Tr(avr)
    K_pa = PSY.get_K_pa(avr) #k_pa>0
    K_ia = PSY.get_K_ia(avr)
    K_da = PSY.get_K_da(avr)
    T_da = PSY.get_T_da(avr)
    Va_min, Va_max = PSY.get_Va_lim(avr)
    K_ff = PSY.get_K_ff(avr)
    K_m = PSY.get_K_m(avr)
    K_ci = PSY.get_K_ci(avr) #K_cl in pss
    K_lr = PSY.get_K_lr(avr)
    I_lr = PSY.get_I_lr(avr)
    Vr_min, Vr_max = PSY.get_Vr_lim(avr)
    Kg = PSY.get_Kg(avr)
    Tg = PSY.get_Tg(avr) #T_g>0

    #Compute block derivatives
    _, dVm_dt = low_pass_mass_matrix(V_th, Vm, 1.0, Tr)
    pid_input = V_ref + Vs - Vm
    pi_out, dx_i = pi_block_nonwindup(pid_input, x_i, K_pa, K_ia, Va_min, Va_max)
    pd_out, dx_d = high_pass_mass_matrix(pid_input, x_d, K_da, T_da)
    Va = pi_out + pd_out

    ff_out = ((Va - Vg) * K_m) + (K_ff * Va)
    V_r1 = max(((I_lr * K_ci) - Xad_Ifd) * K_lr, Vr_min)
    V_r2 = clamp(ff_out, Vr_min, Vr_max)
    V_r = min(V_r1, V_r2)
    E_fd = V_r * Vm
    _, dVg = low_pass(E_fd, Vg, Kg, Tg)

    #Compute 4 States AVR ODE:
    output_ode[local_ix[1]] = dVm_dt
    output_ode[local_ix[2]] = dx_i
    output_ode[local_ix[3]] = dx_d
    output_ode[local_ix[4]] = dVg

    #Update inner_vars
    inner_vars[Vf_var] = E_fd
    return
end

function mdl_avr_ode!(
    device_states::AbstractArray,
    output_ode::AbstractArray,
    inner_vars::AbstractArray,
    dynamic_device::DynamicWrapper{PSY.DynamicGenerator{M, S, PSY.ST8C, TG, P}},
    h,
    t,
) where {M <: PSY.Machine, S <: PSY.Shaft, TG <: PSY.TurbineGov, P <: PSY.PSS}

    #Obtain references
    V_ref = get_V_ref(dynamic_device)

    #Obtain indices for component w/r to device
    local_ix = get_local_state_ix(dynamic_device, PSY.ST8C)

    #Define inner states for component
    internal_states = @view device_states[local_ix]
    Vm = internal_states[1] # Sensed Voltage
    x_a1 = internal_states[2] # Regulator Integrator
    x_a2 = internal_states[3] # Field Regulator
    x_a3 = internal_states[4] # Controller Integrator
    x_a4 = internal_states[5] # Regulator Feedback

    #Define external states for device
    Vt = sqrt(inner_vars[VR_gen_var]^2 + inner_vars[VI_gen_var]^2) # machine's terminal voltage
    Vs = inner_vars[V_pss_var] # PSS output 
    Ifd = inner_vars[Xad_Ifd_var] # machine's field current in exciter base 

    #Get parameters
    avr = PSY.get_avr(dynamic_device)
    SW1_Flag = PSY.get_SW1_flags(avr)
    if SW1_Flag == 1
        error("Source from generator terminal voltage not supported.")
    end
    Tr = PSY.get_Tr(avr)
    K_pr = PSY.get_K_pr(avr)
    K_ir = PSY.get_K_ir(avr)
    Vpi_min, Vpi_max = PSY.get_Vpi_lim(avr)
    K_pa = PSY.get_K_pa(avr)
    K_ia = PSY.get_K_ia(avr)
    Va_min, Va_max = PSY.get_Va_lim(avr)
    K_a = PSY.get_K_a(avr)
    T_a = PSY.get_T_a(avr)
    Vr_min, Vr_max = PSY.get_Vr_lim(avr)
    K_f = PSY.get_K_f(avr)
    T_f = PSY.get_T_f(avr)
    K_c1 = PSY.get_K_c1(avr)
    K_p = PSY.get_K_p(avr)
    K_i2 = PSY.get_K_i2(avr)
    VB1_max = PSY.get_VB1_max(avr)

    if K_i2 != 0.0
        error("Feedforward Current for AVR ST8C not implemented yet.")
    end
    #TODO: Implement Terminal Current FF for AVR
    V_b2 = 0.0
    #TODO: Implement Voltage Compensation if needed
    V_e = K_p

    # Compute block derivatives
    _, dVm_dt = low_pass_mass_matrix(Vt, Vm, 1.0, Tr)
    V_pi_in = V_ref + Vs - Vm
    Ifd_ref, dxa1_dt = pi_block_nonwindup(V_pi_in, x_a1, K_pr, K_ir, Vpi_min, Vpi_max)
    Ifd_diff = Ifd_ref - x_a4
    pi_out, dxa2_dt = pi_block_nonwindup(Ifd_diff, x_a2, K_pa, K_ia, Va_min, Va_max)
    _, dxa3_dt = low_pass_nonwindup_mass_matrix(pi_out, x_a3, K_a, T_a, Vr_min, Vr_max)
    _, dxa4_dt = low_pass_mass_matrix(Ifd, x_a4, K_f, T_f)

    # Compute V_b1
    I_n1 = K_c1 * Ifd / V_e
    F_ex = rectifier_function(I_n1)
    V_b1 = min(F_ex * V_e, VB1_max)

    #Compute 5 States AVR ODE:
    output_ode[local_ix[1]] = dVm_dt
    output_ode[local_ix[2]] = dxa1_dt
    output_ode[local_ix[3]] = dxa2_dt
    output_ode[local_ix[4]] = dxa3_dt
    output_ode[local_ix[5]] = dxa4_dt

    #Update inner_vars
    Efd = V_b1 * x_a3 + V_b2
    inner_vars[Vf_var] = Efd
    return
end
