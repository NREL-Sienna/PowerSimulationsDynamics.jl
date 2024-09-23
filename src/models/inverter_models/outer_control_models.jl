function mass_matrix_outer_entries!(
    mass_matrix,
    outer_control::O,
    global_index::Base.ImmutableDict{Symbol, Int64},
) where {O <: PSY.OuterControl}
    @debug "Using default mass matrix entries $O"
end

#####################################################
### Auxiliary ODE calculations via Flags dispatch ###
#####################################################

### Active Controllers ###

#Freq_Flag = 1
function _mdl_ode_RE_active_controller_AB!(
    active_controller_ode::AbstractArray{<:ACCEPTED_REAL_TYPES},
    active_controller_states::AbstractArray{<:ACCEPTED_REAL_TYPES},
    p::AbstractArray{<:ACCEPTED_REAL_TYPES},
    p_elec_out::ACCEPTED_REAL_TYPES,
    ω_sys::ACCEPTED_REAL_TYPES,
    Vt_filt::ACCEPTED_REAL_TYPES,
    ::Val{1},
    active_power_control::PSY.ActiveRenewableControllerAB,
    dynamic_device::DynamicWrapper{
        PSY.DynamicInverter{
            C,
            PSY.OuterControl{
                PSY.ActiveRenewableControllerAB,
                PSY.ReactiveRenewableControllerAB,
            },
            IC,
            DC,
            P,
            F,
            L,
        },
    },
    inner_vars::AbstractVector,
) where {
    C <: PSY.Converter,
    IC <: PSY.InnerControl,
    DC <: PSY.DCSource,
    P <: PSY.FrequencyEstimator,
    F <: PSY.Filter,
    L <: Union{Nothing, PSY.OutputCurrentLimiter},
}

    #Obtain external parameters
    p_ref = p[:refs][:P_ref]
    ω_ref = p[:refs][:ω_ref]
    # To do: Obtain proper frequency for a plant. For now using the system frequency.
    ω_plant = ω_sys

    #Obtain additional Active Power Controller parameters   
    params = p[:params][:OuterControl][:ActivePowerControl]
    K_pg = params[:K_pg]
    K_ig = params[:K_ig]
    T_p = params[:T_p]
    fe_min = params[:fe_lim][:min]
    fe_max = params[:fe_lim][:max]
    P_min = params[:P_lim][:min]
    P_max = params[:P_lim][:max]
    T_g = params[:T_g]
    D_dn = params[:D_dn]
    D_up = params[:D_up]
    dP_min = params[:dP_lim][:min]
    dP_max = params[:dP_lim][:max]
    P_min_inner = params[:P_lim_inner][:min]
    P_max_inner = params[:P_lim_inner][:max]
    T_pord = params[:T_pord]

    # Not considered parameters because not named tuple
    fdbd1, fdbd2 = PSY.get_fdbd_pnts(active_power_control)

    #Define internal states for outer control
    p_flt = active_controller_states[1]
    ξ_P = active_controller_states[2]
    p_ext = active_controller_states[3]
    p_ord = active_controller_states[4]

    #Compute additional terms
    f_err = deadband_function(ω_ref - ω_plant, fdbd1, fdbd2)
    p_droop = min(D_dn * f_err, 0.0) + max(D_up * f_err, 0.0)
    p_err = clamp(p_ref + p_droop - p_flt, fe_min, fe_max)

    #Compute block derivatives REPCA
    _, dpflt_dt = low_pass(p_elec_out, p_flt, 1.0, T_p)
    P_pi, dξP_dt = pi_block_nonwindup(p_err, ξ_P, K_pg, K_ig, P_min, P_max)
    _, dpext_dt = low_pass(P_pi, p_ext, 1.0, T_g)

    #Compute block derivatives REECB
    p_ord_sat, dpord_dt = low_pass_nonwindup_ramp_limits(
        p_ext,
        p_ord,
        1.0,
        T_pord,
        P_min_inner,
        P_max_inner,
        dP_min,
        dP_max,
    )

    #Update ODEs
    active_controller_ode[1] = dpflt_dt
    active_controller_ode[2] = dξP_dt
    active_controller_ode[3] = dpext_dt
    active_controller_ode[4] = dpord_dt

    #Update Inner Vars: Ioc_pcmd
    inner_vars[Id_oc_var] = p_ord_sat / max(Vt_filt, VOLTAGE_DIVISION_LOWER_BOUND)
    return
end

#Freq_Flag = 0
function _mdl_ode_RE_active_controller_AB!(
    active_controller_ode::AbstractArray{<:ACCEPTED_REAL_TYPES},
    active_controller_states::AbstractArray{<:ACCEPTED_REAL_TYPES},
    p::AbstractArray{<:ACCEPTED_REAL_TYPES},
    p_elec_out::ACCEPTED_REAL_TYPES,
    ω_sys::ACCEPTED_REAL_TYPES,
    Vt_filt::ACCEPTED_REAL_TYPES,
    ::Val{0},
    active_power_control::PSY.ActiveRenewableControllerAB,
    dynamic_device::DynamicWrapper{
        PSY.DynamicInverter{
            C,
            PSY.OuterControl{
                PSY.ActiveRenewableControllerAB,
                PSY.ReactiveRenewableControllerAB,
            },
            IC,
            DC,
            P,
            F,
            L,
        },
    },
    inner_vars::AbstractVector,
) where {
    C <: PSY.Converter,
    IC <: PSY.InnerControl,
    DC <: PSY.DCSource,
    P <: PSY.FrequencyEstimator,
    F <: PSY.Filter,
    L <: Union{Nothing, PSY.OutputCurrentLimiter},
}

    #Obtain external parameters
    p_ref = p[:refs][:P_ref]
    #Obtain additional Active Power Controller parameters
    T_pord = p[:params][:OuterControl][:ActivePowerControl][:T_pord]

    #Define internal states for outer control
    p_ord = active_controller_states[1]

    #Compute block derivatives
    _, dpord_dt = low_pass(p_ref, p_ord, 1.0, T_pord)

    #Update ODE
    active_controller_ode[1] = dpord_dt

    #Update Inner Vars: Ioc_pcmd
    inner_vars[Id_oc_var] = p_ord / max(Vt_filt, VOLTAGE_DIVISION_LOWER_BOUND)
    return
end

### Reactive Controllers ###

# Order of Flags are:
# Ref_Flag, PF_Flag, V_Flag

# VC_Flag == N/A && Ref_Flag == 0 && PF_Flag == 0 && V_Flag == 1
# Named: Fixed Plant Q: Not compliant if Q_Flag = 0 from Inner Control
# Named: Plant Q and Local Q/V if Q_Flag = 1 (compliant).
function _mdl_ode_RE_reactive_controller_AB!(
    reactive_controller_ode::AbstractArray{<:ACCEPTED_REAL_TYPES},
    reactive_controller_states::AbstractArray{<:ACCEPTED_REAL_TYPES},
    p::AbstractArray{<:ACCEPTED_REAL_TYPES},
    q_elec_out::ACCEPTED_REAL_TYPES,
    Vt_filt::ACCEPTED_REAL_TYPES,
    ::Val{0},
    ::Val{0},
    ::Val{1},
    reactive_power_control::PSY.ReactiveRenewableControllerAB,
    dynamic_device::DynamicWrapper{
        PSY.DynamicInverter{
            C,
            PSY.OuterControl{
                PSY.ActiveRenewableControllerAB,
                PSY.ReactiveRenewableControllerAB,
            },
            IC,
            DC,
            P,
            F,
            L,
        },
    },
    inner_vars::AbstractVector,
) where {
    C <: PSY.Converter,
    IC <: PSY.InnerControl,
    DC <: PSY.DCSource,
    P <: PSY.FrequencyEstimator,
    F <: PSY.Filter,
    L <: Union{Nothing, PSY.OutputCurrentLimiter},
}

    #Obtain external parameters
    q_ref = p[:refs][:Q_ref]
    params = p[:params][:OuterControl][:ReactivePowerControl]
    # Get Reactive Controller parameters
    T_fltr = params[:T_fltr]
    K_p = params[:K_p]
    K_i = params[:K_i]
    T_ft = params[:T_ft]
    T_fv = params[:T_fv]
    e_min = params[:e_lim][:min]
    e_max = params[:e_lim][:max]
    dbd1 = params[:dbd_pnts1]
    dbd2 = params[:dbd_pnts2]
    Q_min = params[:Q_lim][:min]
    Q_max = params[:Q_lim][:max]
    # V_frz not implemented yet
    Q_min_inner = params[:Q_lim_inner][:min]
    Q_max_inner = params[:Q_lim_inner][:max]
    V_min = params[:V_lim][:min]
    V_max = params[:V_lim][:max]
    K_qp = params[:K_qp]
    K_qi = params[:K_qi]

    #Define internal states for Reactive Control
    q_flt = reactive_controller_states[1]
    ξq_oc = reactive_controller_states[2]
    q_LL = reactive_controller_states[3]
    ξ_Q = reactive_controller_states[4]

    #Compute additional variables
    q_err = clamp(deadband_function(q_ref - q_flt, dbd1, dbd2), e_min, e_max)

    #Compute block derivatives REPCA
    _, dqflt_dt = low_pass(q_elec_out, q_flt, 1.0, T_fltr)
    Q_pi_sat, dξqoc_dt = pi_block_nonwindup(q_err, ξq_oc, K_p, K_i, Q_min, Q_max)
    Q_ext, dqLL_dt = lead_lag(Q_pi_sat, q_LL, 1.0, T_ft, T_fv)

    #Compute block derivatives REECB
    V_pi_in = clamp(Q_ext, Q_min_inner, Q_max_inner) - q_elec_out
    V_pi_sat, dξQ_dt = pi_block_nonwindup(V_pi_in, ξ_Q, K_qp, K_qi, V_min, V_max)

    #Update ODEs
    reactive_controller_ode[1] = dqflt_dt
    reactive_controller_ode[2] = dξqoc_dt
    reactive_controller_ode[3] = dqLL_dt
    reactive_controller_ode[4] = dξQ_dt

    #Update Inner Vars
    inner_vars[V_oc_var] = V_pi_sat - Vt_filt
    inner_vars[Iq_oc_var] = Q_ext / max(Vt_filt, VOLTAGE_DIVISION_LOWER_BOUND)
    return
end

# VC_Flag == N/A && Ref_Flag == 0 && PF_Flag == 0 && V_Flag == 0
# Fixed Plant Q if Q_Flag = 0 from Inner Control (not Compliant from CAISO requirements)
function _mdl_ode_RE_reactive_controller_AB!(
    reactive_controller_ode::AbstractArray{<:ACCEPTED_REAL_TYPES},
    reactive_controller_states::AbstractArray{<:ACCEPTED_REAL_TYPES},
    p::AbstractArray{<:ACCEPTED_REAL_TYPES},
    q_elec_out::ACCEPTED_REAL_TYPES,
    Vt_filt::ACCEPTED_REAL_TYPES,
    ::Val{0},
    ::Val{0},
    ::Val{0},
    reactive_power_control::PSY.ReactiveRenewableControllerAB,
    dynamic_device::DynamicWrapper{
        PSY.DynamicInverter{
            C,
            PSY.OuterControl{
                PSY.ActiveRenewableControllerAB,
                PSY.ReactiveRenewableControllerAB,
            },
            IC,
            DC,
            P,
            F,
            L,
        },
    },
    inner_vars::AbstractVector,
) where {
    C <: PSY.Converter,
    IC <: PSY.InnerControl,
    DC <: PSY.DCSource,
    P <: PSY.FrequencyEstimator,
    F <: PSY.Filter,
    L <: Union{Nothing, PSY.OutputCurrentLimiter},
}
    #Obtain external parameters
    q_ref = p[:refs][:Q_ref]
    # Get Reactive Controller parameters
    params = p[:params][:OuterControl][:ReactivePowerControl]
    T_fltr = params[:T_fltr]
    K_p = params[:K_p]
    K_i = params[:K_i]
    T_ft = params[:T_ft]
    T_fv = params[:T_fv]
    e_min = params[:e_lim][:min]
    e_max = params[:e_lim][:max]
    dbd1 = params[:dbd_pnts1]
    dbd2 = params[:dbd_pnts2]
    Q_min = params[:Q_lim][:min]
    Q_max = params[:Q_lim][:max]
    # V_frz not implemented yet
    #Define internal states for Reactive Control
    q_flt = reactive_controller_states[1]
    ξq_oc = reactive_controller_states[2]
    q_LL = reactive_controller_states[3]

    #Compute additional variables
    q_err = clamp(deadband_function(q_ref - q_flt, dbd1, dbd2), e_min, e_max)

    #Compute block derivatives REPCA
    _, dqflt_dt = low_pass(q_elec_out, q_flt, 1.0, T_fltr)
    Q_pi_sat, dξqoc_dt = pi_block_nonwindup(q_err, ξq_oc, K_p, K_i, Q_min, Q_max)
    Q_ext, dqLL_dt = lead_lag(Q_pi_sat, q_LL, 1.0, T_ft, T_fv)

    #Update ODEs
    reactive_controller_ode[1] = dqflt_dt
    reactive_controller_ode[2] = dξqoc_dt
    reactive_controller_ode[3] = dqLL_dt

    #Update Inner Vars
    inner_vars[V_oc_var] = Q_ext - Vt_filt
    inner_vars[Iq_oc_var] = Q_ext / max(Vt_filt, VOLTAGE_DIVISION_LOWER_BOUND)
    return
end

# VC_Flag == 0 or 1 && Ref_Flag == 1 && PF_Flag == 0 && V_Flag == 1
# Plant V if Q_Flag = 0 (compliant)
# Plant V and Local Q/V if Q_Flag = 1 (compliant)
# Ref_Flag = 1 is used for plant-level voltage control, typically some point of the entire plant
# This depends on where V_reg is being measured. For our purposes, to regulate V_t, the voltage bus that
# the plant is connected, you can do it by setting up K_c = 0 (for VC_Flag = 0) or Rc = Xc = 0 (for VC_Flag = 1)
function _mdl_ode_RE_reactive_controller_AB!(
    reactive_controller_ode,
    reactive_controller_states,
    p,
    q_elec_out,
    Vt_filt,
    ::Val{1},
    ::Val{0},
    ::Val{1},
    reactive_power_control::PSY.ReactiveRenewableControllerAB,
    dynamic_device::DynamicWrapper{
        PSY.DynamicInverter{
            C,
            PSY.OuterControl{
                PSY.ActiveRenewableControllerAB,
                PSY.ReactiveRenewableControllerAB,
            },
            IC,
            DC,
            P,
            F,
            L,
        },
    },
    inner_vars::AbstractVector,
) where {
    C <: PSY.Converter,
    IC <: PSY.InnerControl,
    DC <: PSY.DCSource,
    P <: PSY.FrequencyEstimator,
    F <: PSY.Filter,
    L <: Union{Nothing, PSY.OutputCurrentLimiter},
}
    #Obtain external parameters
    V_ref = p[:refs][:V_ref]

    #Obtain regulated voltage (assumed to be terminal voltage)
    V_reg = sqrt(inner_vars[Vr_inv_var]^2 + inner_vars[Vi_inv_var]^2)
    I_R = inner_vars[Ir_inv_var]
    I_I = inner_vars[Ii_inv_var]

    VC_Flag = PSY.get_VC_Flag(reactive_power_control)

    # Get Reactive Controller parameters
    params = p[:params][:OuterControl][:ReactivePowerControl]
    T_fltr = params[:T_fltr]
    K_p = params[:K_p]
    K_i = params[:K_i]
    T_ft = params[:T_ft]
    T_fv = params[:T_fv]
    e_min = params[:e_lim][:min]
    e_max = params[:e_lim][:max]
    dbd1 = params[:dbd_pnts1]
    dbd2 = params[:dbd_pnts2]
    Q_min = params[:Q_lim][:min]
    Q_max = params[:Q_lim][:max]
    # V_frz not implemented yet
    R_c = params[:R_c]
    X_c = params[:X_c]
    K_c = params[:K_c]
    Q_min_inner = params[:Q_lim_inner][:min]
    Q_max_inner = params[:Q_lim_inner][:max]
    V_min = params[:V_lim][:min]
    V_max = params[:V_lim][:max]
    K_qp = params[:K_qp]
    K_qi = params[:K_qi]

    #Define internal states for Reactive Control
    V_cflt = reactive_controller_states[1]
    ξq_oc = reactive_controller_states[2]
    q_LL = reactive_controller_states[3]
    ξ_Q = reactive_controller_states[4]

    #Compute input to the compensated voltage filter
    if VC_Flag == 0
        V_flt_input = V_reg + K_c * q_elec_out
    else
        # Calculate compensated voltage: | V_reg - (R_c + jX_c)(I_r + jI_i) |
        V_flt_input = sqrt(
            V_reg^2 +
            2 * V_reg * (I_I * X_c - I_R * R_c) +
            (I_I^2 + I_R^2) * (R_c^2 + X_c^2),
        )
    end
    #Q error - PI controller
    q_err = clamp(deadband_function(V_ref - V_cflt, dbd1, dbd2), e_min, e_max)

    #Compute block derivatives REPCA
    _, dVcflt_dt = low_pass(V_flt_input, V_cflt, 1.0, T_fltr)
    Q_pi_sat, dξqoc_dt = pi_block_nonwindup(q_err, ξq_oc, K_p, K_i, Q_min, Q_max)
    Q_ext, dqLL_dt = lead_lag(Q_pi_sat, q_LL, 1.0, T_ft, T_fv)

    #Compute block derivatives REECB
    V_pi_in = clamp(Q_ext, Q_min_inner, Q_max_inner) - q_elec_out
    V_pi_sat, dξQ_dt = pi_block_nonwindup(V_pi_in, ξ_Q, K_qp, K_qi, V_min, V_max)

    #Update ODEs
    reactive_controller_ode[1] = dVcflt_dt
    reactive_controller_ode[2] = dξqoc_dt
    reactive_controller_ode[3] = dqLL_dt
    reactive_controller_ode[4] = dξQ_dt

    #Update Inner Vars
    inner_vars[V_oc_var] = V_pi_sat - Vt_filt
    inner_vars[Iq_oc_var] = Q_ext / max(Vt_filt, VOLTAGE_DIVISION_LOWER_BOUND)
    return
end

# VC_Flag == 0 or 1 && Ref_Flag == 1 && PF_Flag == 0 && V_Flag == 0
# Plant V if Q_Flag = 0 (compliant)
# Plant V and Local V if Q_Flag = 1 (compliant)
# Ref_Flag = 1 is used for plant-level voltage control, typically some point of the entire plant
# This depends on where V_reg is being measured. For our purposes, to regulate V_t, the voltage bus that
# the plant is connected, you can do it by setting up K_c = 0 (for VC_Flag = 0) or Rc = Xc = 0 (for VC_Flag = 1)
# CURRENTLY NOT WORKING with Q_Flag = 1
function _mdl_ode_RE_reactive_controller_AB!(
    reactive_controller_ode,
    reactive_controller_states,
    p,
    q_elec_out,
    Vt_filt,
    ::Val{1},
    ::Val{0},
    ::Val{0},
    reactive_power_control::PSY.ReactiveRenewableControllerAB,
    dynamic_device::DynamicWrapper{
        PSY.DynamicInverter{
            C,
            PSY.OuterControl{
                PSY.ActiveRenewableControllerAB,
                PSY.ReactiveRenewableControllerAB,
            },
            IC,
            DC,
            P,
            F,
            L,
        },
    },
    inner_vars::AbstractVector,
) where {
    C <: PSY.Converter,
    IC <: PSY.InnerControl,
    DC <: PSY.DCSource,
    P <: PSY.FrequencyEstimator,
    F <: PSY.Filter,
    L <: Union{Nothing, PSY.OutputCurrentLimiter},
}
    #Obtain external parameters
    V_ref = p[:refs][:V_ref]

    #Obtain regulated voltage (assumed to be terminal voltage)
    V_reg = sqrt(inner_vars[Vr_inv_var]^2 + inner_vars[Vi_inv_var]^2)
    I_R = inner_vars[Ir_inv_var]
    I_I = inner_vars[Ii_inv_var]

    # Get Reactive Controller parameters
    VC_Flag = PSY.get_VC_Flag(reactive_power_control)

    params = p[:params][:OuterControl][:ReactivePowerControl]
    T_fltr = params[:T_fltr]
    K_p = params[:K_p]
    K_i = params[:K_i]
    T_ft = params[:T_ft]
    T_fv = params[:T_fv]
    e_min = params[:e_lim][:min]
    e_max = params[:e_lim][:max]
    dbd1 = params[:dbd_pnts1]
    dbd2 = params[:dbd_pnts2]
    Q_min = params[:Q_lim][:min]
    Q_max = params[:Q_lim][:max]
    # V_frz not implemented yet
    R_c = params[:R_c]
    X_c = params[:X_c]
    K_c = params[:K_c]
    V_min = params[:V_lim][:min]
    V_max = params[:V_lim][:max]

    #Define internal states for Reactive Control
    V_cflt = reactive_controller_states[1]
    ξq_oc = reactive_controller_states[2]
    q_LL = reactive_controller_states[3]

    #Compute input to the compensated voltage filter
    if VC_Flag == 0
        V_flt_input = V_reg + K_c * q_elec_out
    else
        # Calculate compensated voltage: | V_reg - (R_c + jX_c)(I_r + jI_i) |
        V_flt_input = sqrt(
            V_reg^2 +
            2 * V_reg * (I_I * X_c - I_R * R_c) +
            (I_I^2 + I_R^2) * (R_c^2 + X_c^2),
        )
    end
    #Q error - PI controller
    q_err = clamp(deadband_function(V_ref - V_cflt, dbd1, dbd2), e_min, e_max)

    #Compute block derivatives REPCA
    _, dVcflt_dt = low_pass(V_flt_input, V_cflt, 1.0, T_fltr)
    Q_pi_sat, dξqoc_dt = pi_block_nonwindup(q_err, ξq_oc, K_p, K_i, Q_min, Q_max)
    Q_ext, dqLL_dt = lead_lag(Q_pi_sat, q_LL, 1.0, T_ft, T_fv)

    #Skip PI voltage inner block
    V_pi_sat = clamp(Q_ext, V_min, V_max)

    #Update ODEs
    reactive_controller_ode[1] = dVcflt_dt
    reactive_controller_ode[2] = dξqoc_dt
    reactive_controller_ode[3] = dqLL_dt

    #Update Inner Vars
    inner_vars[V_oc_var] = V_pi_sat - Vt_filt
    inner_vars[Iq_oc_var] = Q_ext / max(Vt_filt, VOLTAGE_DIVISION_LOWER_BOUND)
    return
end

############################################
### ODE calculations via device dispatch ###
############################################

function mdl_outer_ode!(
    device_states::AbstractArray{<:ACCEPTED_REAL_TYPES},
    output_ode::AbstractArray{<:ACCEPTED_REAL_TYPES},
    p::AbstractArray{<:ACCEPTED_REAL_TYPES},
    inner_vars::AbstractArray{<:ACCEPTED_REAL_TYPES},
    ω_sys::ACCEPTED_REAL_TYPES,
    dynamic_device::DynamicWrapper{
        PSY.DynamicInverter{
            C,
            PSY.OuterControl{PSY.VirtualInertia, PSY.ReactivePowerDroop},
            IC,
            DC,
            P,
            F,
            L,
        },
    },
    h,
    t,
) where {
    C <: PSY.Converter,
    IC <: PSY.InnerControl,
    DC <: PSY.DCSource,
    P <: PSY.FrequencyEstimator,
    F <: PSY.Filter,
    L <: Union{Nothing, PSY.OutputCurrentLimiter},
}

    #Obtain external states inputs for component
    external_ix = get_input_port_ix(
        dynamic_device,
        PSY.OuterControl{PSY.VirtualInertia, PSY.ReactivePowerDroop},
    )
    Vr_filter = device_states[external_ix[1]]
    Vi_filter = device_states[external_ix[2]]
    Ir_filter = device_states[external_ix[3]]
    Ii_filter = device_states[external_ix[4]]

    #Obtain inner variables for component
    ω_pll = inner_vars[ω_freq_estimator_var]

    Ta = p[:params][:OuterControl][:ActivePowerControl][:Ta]
    kd = p[:params][:OuterControl][:ActivePowerControl][:kd]
    kω = p[:params][:OuterControl][:ActivePowerControl][:kω]
    kq = p[:params][:OuterControl][:ReactivePowerControl][:kq]
    ωf = p[:params][:OuterControl][:ReactivePowerControl][:ωf]

    q_ref = p[:refs][:Q_ref]
    V_ref = p[:refs][:V_ref]
    ω_ref = p[:refs][:ω_ref]
    p_ref = p[:refs][:P_ref]

    f0 = get_system_base_frequency(dynamic_device)
    ωb = 2 * pi * f0 #Rated angular frequency

    #Obtain indices for component w/r to device
    local_ix = get_local_state_ix(
        dynamic_device,
        PSY.OuterControl{PSY.VirtualInertia, PSY.ReactivePowerDroop},
    )

    #Define internal states for Outer Control
    internal_states = @view device_states[local_ix]
    θ_oc = internal_states[1]
    ω_oc = internal_states[2]
    qm = internal_states[3]

    #Obtain additional expressions
    p_elec_out = Ir_filter * Vr_filter + Ii_filter * Vi_filter
    q_elec_out = -Ii_filter * Vr_filter + Ir_filter * Vi_filter

    #Compute block derivatives
    _, dqm_dt = low_pass(q_elec_out, qm, 1.0, 1.0 / ωf)

    #Compute 3 states ODEs
    output_ode[local_ix[1]] = ωb * (ω_oc - ω_sys)
    output_ode[local_ix[2]] =
        (p_ref / Ta - p_elec_out / Ta - kd * (ω_oc - ω_pll) / Ta - kω * (ω_oc - ω_ref) / Ta)
    output_ode[local_ix[3]] = dqm_dt

    #Update inner vars
    inner_vars[θ_oc_var] = θ_oc
    inner_vars[ω_oc_var] = ω_oc
    inner_vars[V_oc_var] = V_ref + kq * (q_ref - qm)
    return
end

function mdl_outer_ode!(
    device_states::AbstractArray{<:ACCEPTED_REAL_TYPES},
    output_ode::AbstractArray{<:ACCEPTED_REAL_TYPES},
    p::AbstractArray{<:ACCEPTED_REAL_TYPES},
    inner_vars::AbstractArray{<:ACCEPTED_REAL_TYPES},
    ω_sys::ACCEPTED_REAL_TYPES,
    dynamic_device::DynamicWrapper{
        PSY.DynamicInverter{
            C,
            PSY.OuterControl{PSY.ActivePowerDroop, PSY.ReactivePowerDroop},
            IC,
            DC,
            P,
            F,
            L,
        },
    },
    h,
    t,
) where {
    C <: PSY.Converter,
    IC <: PSY.InnerControl,
    DC <: PSY.DCSource,
    P <: PSY.FrequencyEstimator,
    F <: PSY.Filter,
    L <: Union{Nothing, PSY.OutputCurrentLimiter},
}

    #Obtain external states inputs for component
    external_ix = get_input_port_ix(
        dynamic_device,
        PSY.OuterControl{PSY.ActivePowerDroop, PSY.ReactivePowerDroop},
    )
    Vr_filter = device_states[external_ix[1]]
    Vi_filter = device_states[external_ix[2]]
    Ir_filter = device_states[external_ix[3]]
    Ii_filter = device_states[external_ix[4]]

    #Get Active Power Controller parameters
    Rp = p[:params][:OuterControl][:ActivePowerControl][:Rp]   #Droop Gain
    ωz = p[:params][:OuterControl][:ActivePowerControl][:ωz]  #Frequency cutoff frequency

    f0 = get_system_base_frequency(dynamic_device)
    ωb = 2 * pi * f0 #Rated angular frequency

    #Get Reactive Power Controller parameters
    kq = p[:params][:OuterControl][:ReactivePowerControl][:kq]
    ωf = p[:params][:OuterControl][:ReactivePowerControl][:ωf]

    #Obtain external parameters
    p_ref = p[:refs][:P_ref]
    ω_ref = p[:refs][:ω_ref]
    V_ref = p[:refs][:V_ref]
    q_ref = p[:refs][:Q_ref]

    #Obtain indices for component w/r to device
    local_ix = get_local_state_ix(
        dynamic_device,
        PSY.OuterControl{PSY.ActivePowerDroop, PSY.ReactivePowerDroop},
    )

    #Define internal states for frequency estimator
    internal_states = @view device_states[local_ix]
    θ_oc = internal_states[1]
    pm = internal_states[2]
    qm = internal_states[3]

    #Obtain additional expressions
    p_elec_out = Ir_filter * Vr_filter + Ii_filter * Vi_filter
    q_elec_out = -Ii_filter * Vr_filter + Ir_filter * Vi_filter

    #Compute Frequency from Droop
    ω_oc = ω_ref + Rp * (p_ref - pm)

    #Compute block derivatives
    _, dpm_dt = low_pass(p_elec_out, pm, 1.0, 1.0 / ωz)
    _, dqm_dt = low_pass(q_elec_out, qm, 1.0, 1.0 / ωf)

    outer_control = PSY.get_outer_control(dynamic_device)
    ext = PSY.get_ext(outer_control)
    bool_val = get(ext, "is_not_reference", 1.0)

    #Compute 3 states ODEs
    output_ode[local_ix[1]] = bool_val * ωb * (ω_oc - ω_sys)
    output_ode[local_ix[2]] = dpm_dt
    output_ode[local_ix[3]] = dqm_dt

    #Update inner vars
    inner_vars[θ_oc_var] = θ_oc
    inner_vars[ω_oc_var] = ω_oc
    inner_vars[V_oc_var] = V_ref + kq * (q_ref - qm)
    return
end

function mdl_outer_ode!(
    device_states::AbstractArray{<:ACCEPTED_REAL_TYPES},
    output_ode::AbstractArray{<:ACCEPTED_REAL_TYPES},
    p::AbstractArray{<:ACCEPTED_REAL_TYPES},
    inner_vars::AbstractArray{<:ACCEPTED_REAL_TYPES},
    ω_sys::ACCEPTED_REAL_TYPES,
    dynamic_device::DynamicWrapper{
        PSY.DynamicInverter{
            C,
            PSY.OuterControl{PSY.ActiveVirtualOscillator, PSY.ReactiveVirtualOscillator},
            IC,
            DC,
            P,
            F,
            L,
        },
    },
    h,
    t,
) where {
    C <: PSY.Converter,
    IC <: PSY.InnerControl,
    DC <: PSY.DCSource,
    P <: PSY.FrequencyEstimator,
    F <: PSY.Filter,
    L <: Union{Nothing, PSY.OutputCurrentLimiter},
}

    #Obtain external states inputs for component
    external_ix = get_input_port_ix(
        dynamic_device,
        PSY.OuterControl{PSY.ActiveVirtualOscillator, PSY.ReactiveVirtualOscillator},
    )
    Vr_filter = device_states[external_ix[1]]
    Vi_filter = device_states[external_ix[2]]
    Ir_filter = device_states[external_ix[3]]
    Ii_filter = device_states[external_ix[4]]

    #Get Active Power Controller parameters
    params_active = p[:params][:OuterControl][:ActivePowerControl]
    k1 = params_active[:k1]
    ψ = params_active[:ψ]
    f0 = get_system_base_frequency(dynamic_device)
    ωb = 2 * pi * f0 #Rated angular frequency

    #Get Reactive Power Controller parameters
    params_reactive = p[:params][:OuterControl][:ReactivePowerControl]
    k2 = params_reactive[:k2]

    #Obtain external parameters
    p_ref = p[:refs][:P_ref]
    V_ref = p[:refs][:V_ref]
    q_ref = p[:refs][:Q_ref]

    #Obtain indices for component w/r to device
    local_ix = get_local_state_ix(
        dynamic_device,
        PSY.OuterControl{PSY.ActiveVirtualOscillator, PSY.ReactiveVirtualOscillator},
    )

    #Define internal states for frequency estimator
    internal_states = @view device_states[local_ix]
    θ_oc = internal_states[1]
    E_oc = internal_states[2]

    #Obtain additional expressions
    p_elec_out = Ir_filter * Vr_filter + Ii_filter * Vi_filter
    q_elec_out = -Ii_filter * Vr_filter + Ir_filter * Vi_filter

    #Compute Frequency from VOC
    γ = ψ - pi / 2
    ω_oc =
        ω_sys +
        (k1 / E_oc^2) * (cos(γ) * (p_ref - p_elec_out) + sin(γ) * (q_ref - q_elec_out))

    #Compute voltage derivatives
    dEoc_dt =
        ωb * (
            (k1 / E_oc) * (-sin(γ) * (p_ref - p_elec_out) + cos(γ) * (q_ref - q_elec_out)) +
            k2 * (V_ref^2 - E_oc^2) * E_oc
        )
    outer_control = PSY.get_outer_control(dynamic_device)
    ext = PSY.get_ext(outer_control)
    bool_val = get(ext, "is_not_reference", 1.0)

    #Compute 2 states ODEs
    output_ode[local_ix[1]] = bool_val * ωb * (ω_oc - ω_sys)
    output_ode[local_ix[2]] = dEoc_dt

    #Update inner vars
    inner_vars[θ_oc_var] = θ_oc
    inner_vars[ω_oc_var] = ω_oc
    inner_vars[V_oc_var] = E_oc
    return
end

function mdl_outer_ode!(
    device_states::AbstractArray{<:ACCEPTED_REAL_TYPES},
    output_ode::AbstractArray{<:ACCEPTED_REAL_TYPES},
    p::AbstractArray{<:ACCEPTED_REAL_TYPES},
    inner_vars::AbstractArray{<:ACCEPTED_REAL_TYPES},
    ::ACCEPTED_REAL_TYPES,
    dynamic_device::DynamicWrapper{
        PSY.DynamicInverter{
            C,
            PSY.OuterControl{PSY.ActivePowerPI, PSY.ReactivePowerPI},
            IC,
            DC,
            P,
            F,
            L,
        },
    },
    h,
    t,
) where {
    C <: PSY.Converter,
    IC <: PSY.InnerControl,
    DC <: PSY.DCSource,
    P <: PSY.FrequencyEstimator,
    F <: PSY.Filter,
    L <: Union{Nothing, PSY.OutputCurrentLimiter},
}

    #Obtain external states inputs for component
    external_ix = get_input_port_ix(
        dynamic_device,
        PSY.OuterControl{PSY.ActivePowerPI, PSY.ReactivePowerPI},
    )
    Vr_filter = device_states[external_ix[1]]
    Vi_filter = device_states[external_ix[2]]
    Ir_filter = device_states[external_ix[3]]
    Ii_filter = device_states[external_ix[4]]

    #Obtain inner variables for component
    θ_pll = inner_vars[θ_freq_estimator_var]
    ω_pll = inner_vars[ω_freq_estimator_var]

    #Get Active Power Controller parameters
    Kp_p = p[:params][:OuterControl][:ActivePowerControl][:Kp_p]
    Ki_p = p[:params][:OuterControl][:ActivePowerControl][:Ki_p]
    ωz = p[:params][:OuterControl][:ActivePowerControl][:ωz]

    #Get Reactive Power Controller parameters
    Kp_q = p[:params][:OuterControl][:ReactivePowerControl][:Kp_q]
    Ki_q = p[:params][:OuterControl][:ReactivePowerControl][:Ki_q]
    ωf = p[:params][:OuterControl][:ReactivePowerControl][:ωf]

    #Obtain external parameters
    p_ref = p[:refs][:P_ref]
    q_ref = p[:refs][:Q_ref]

    #Obtain indices for component w/r to device
    local_ix = get_local_state_ix(
        dynamic_device,
        PSY.OuterControl{PSY.ActivePowerPI, PSY.ReactivePowerPI},
    )

    #Define internal states for outer control
    internal_states = @view device_states[local_ix]
    σp_oc = internal_states[1]
    p_oc = internal_states[2]
    σq_oc = internal_states[3]
    q_oc = internal_states[4]

    #Obtain additional expressions
    p_elec_out = Ir_filter * Vr_filter + Ii_filter * Vi_filter
    q_elec_out = -Ii_filter * Vr_filter + Ir_filter * Vi_filter

    #Compute block derivatives
    Iq_pi, dσpoc_dt = pi_block(p_ref - p_oc, σp_oc, Kp_p, Ki_p)
    _, dpoc_dt = low_pass(p_elec_out, p_oc, 1.0, 1.0 / ωz)
    Id_pi, dσqoc_dt = pi_block(q_ref - q_oc, σq_oc, Kp_q, Ki_q)
    _, dqoc_dt = low_pass(q_elec_out, q_oc, 1.0, 1.0 / ωf)

    #Compute 4 states ODEs
    output_ode[local_ix[1]] = dσpoc_dt
    output_ode[local_ix[2]] = dpoc_dt
    output_ode[local_ix[3]] = dσqoc_dt
    output_ode[local_ix[4]] = dqoc_dt

    #Update inner vars
    inner_vars[θ_oc_var] = θ_pll
    inner_vars[ω_oc_var] = ω_pll
    inner_vars[Iq_oc_var] = Iq_pi
    inner_vars[Id_oc_var] = Id_pi
    return
end

function mdl_outer_ode!(
    device_states::AbstractArray{<:ACCEPTED_REAL_TYPES},
    output_ode::AbstractArray{<:ACCEPTED_REAL_TYPES},
    p::AbstractArray{<:ACCEPTED_REAL_TYPES},
    inner_vars::AbstractArray{<:ACCEPTED_REAL_TYPES},
    ω_sys::ACCEPTED_REAL_TYPES,
    dynamic_device::DynamicWrapper{
        PSY.DynamicInverter{
            C,
            PSY.OuterControl{
                PSY.ActiveRenewableControllerAB,
                PSY.ReactiveRenewableControllerAB,
            },
            IC,
            DC,
            P,
            F,
            L,
        },
    },
    h,
    t,
) where {
    C <: PSY.Converter,
    IC <: PSY.InnerControl,
    DC <: PSY.DCSource,
    P <: PSY.FrequencyEstimator,
    F <: PSY.Filter,
    L <: Union{Nothing, PSY.OutputCurrentLimiter},
}

    #Obtain external states inputs for component
    external_ix = get_input_port_ix(
        dynamic_device,
        PSY.OuterControl{
            PSY.ActiveRenewableControllerAB,
            PSY.ReactiveRenewableControllerAB,
        },
    )
    inner_ctrl = PSY.get_inner_control(dynamic_device)
    inner_ctrl_ix = get_local_state_ix(dynamic_device, typeof(inner_ctrl))
    inner_ctrl_states = @view device_states[inner_ctrl_ix]
    Vt_filt = inner_ctrl_states[1]

    #Monitoring power from other branch not supported.
    Vr_filter = inner_vars[Vr_filter_var]
    Vi_filter = inner_vars[Vi_filter_var]
    Ir_filter = inner_vars[Ir_filter_var]
    Ii_filter = inner_vars[Ii_filter_var]
    p_elec_out = Ir_filter * Vr_filter + Ii_filter * Vi_filter

    Ir_cnv = inner_vars[Ir_cnv_var]
    Ii_cnv = inner_vars[Ii_cnv_var]
    Ir_cap = Ir_filter - Ir_cnv
    Ii_cap = Ii_filter - Ii_cnv
    Q_cap = -Ii_cap * Vr_filter + Ir_cap * Vi_filter
    q_elec_out = -Ii_filter * Vr_filter + Ir_filter * Vi_filter - Q_cap

    #Get Active Power Controller parameters
    outer_control = PSY.get_outer_control(dynamic_device)
    active_power_control = PSY.get_active_power_control(outer_control)
    #Note: Monitoring power from other branch not supported.
    Freq_Flag = PSY.get_Freq_Flag(active_power_control) #Frequency Flag

    #Get Reactive Power Controller parameters
    reactive_power_control = PSY.get_reactive_power_control(outer_control)
    #Note: Monitoring power from other branch not supported.
    Ref_Flag = PSY.get_Ref_Flag(reactive_power_control)
    PF_Flag = PSY.get_PF_Flag(reactive_power_control)
    V_Flag = PSY.get_V_Flag(reactive_power_control)

    #Obtain indices for component w/r to device
    local_ix = get_local_state_ix(
        dynamic_device,
        PSY.OuterControl{
            PSY.ActiveRenewableControllerAB,
            PSY.ReactiveRenewableControllerAB,
        },
    )

    internal_states = @view device_states[local_ix]
    internal_ode = @view output_ode[local_ix]
    active_n_states = PSY.get_n_states(active_power_control)
    reactive_n_states = PSY.get_n_states(reactive_power_control)
    active_ix_range = 1:active_n_states
    reactive_ix_range = (active_n_states + 1):(active_n_states + reactive_n_states)
    active_states = @view internal_states[active_ix_range]
    reactive_states = @view internal_states[reactive_ix_range]
    active_ode = @view internal_ode[active_ix_range]
    reactive_ode = @view internal_ode[reactive_ix_range]

    #Dispatch active power controller
    _mdl_ode_RE_active_controller_AB!(
        active_ode,
        active_states,
        p,
        p_elec_out,
        ω_sys,
        Vt_filt,
        Val(Freq_Flag),
        active_power_control,
        dynamic_device,
        inner_vars,
    )

    #Dispatch reactive power controller
    _mdl_ode_RE_reactive_controller_AB!(
        reactive_ode,
        reactive_states,
        p,
        q_elec_out,
        Vt_filt,
        Val(Ref_Flag),
        Val(PF_Flag),
        Val(V_Flag),
        reactive_power_control,
        dynamic_device,
        inner_vars,
    )
    return
end
