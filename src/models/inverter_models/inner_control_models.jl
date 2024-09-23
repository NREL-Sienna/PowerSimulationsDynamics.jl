function mass_matrix_inner_entries!(
    mass_matrix,
    inner_control::IC,
    global_index::Base.ImmutableDict{Symbol, Int64},
) where {IC <: PSY.InnerControl}
    @debug "Using default mass matrix entries $IC"
    return
end

function mass_matrix_inner_entries!(
    mass_matrix,
    inner_control::PSY.RECurrentControlB,
    global_index::ImmutableDict{Symbol, Int64},
)
    mass_matrix[global_index[:Vt_filt], global_index[:Vt_filt]] =
        PSY.get_T_rv(inner_control)
    if PSY.get_Q_Flag(inner_control) == 0
        mass_matrix[global_index[:I_icv], global_index[:I_icv]] =
            PSY.get_T_iq(inner_control)
    end
    return
end

#####################################################
### Auxiliary ODE calculations via Flags dispatch ###
#####################################################

### Inner Controllers ###

#Q_Flag = 0
function _mdl_ode_RE_inner_controller_B!(
    inner_controller_ode::AbstractArray{<:ACCEPTED_REAL_TYPES},
    inner_controller_states::AbstractArray{<:ACCEPTED_REAL_TYPES},
    inner_controller_parameters::AbstractArray{<:ACCEPTED_REAL_TYPES},
    ::Val{0},
    inner_control::PSY.RECurrentControlB,
    dynamic_device::DynamicWrapper{
        PSY.DynamicInverter{C, O, PSY.RECurrentControlB, DC, P, F, L},
    },
    inner_vars::AbstractVector,
    h,
    t,
) where {
    C <: PSY.Converter,
    O <: PSY.OuterControl,
    DC <: PSY.DCSource,
    P <: PSY.FrequencyEstimator,
    F <: PSY.Filter,
    L <: Union{Nothing, PSY.OutputCurrentLimiter},
}
    #Obtain inner variables for component
    V_t = sqrt(inner_vars[Vr_filter_var]^2 + inner_vars[Vi_filter_var]^2)
    Ip_oc = inner_vars[Id_oc_var]
    Iq_oc = inner_vars[Iq_oc_var]

    #Get Current Controller parameters
    PQ_Flag = PSY.get_PQ_Flag(inner_control)
    T_rv = inner_controller_parameters[:T_rv]
    dbd1 = inner_controller_parameters[:dbd_pnts1]
    dbd2 = inner_controller_parameters[:dbd_pnts2]
    K_qv = inner_controller_parameters[:K_qv]
    I_ql1 = inner_controller_parameters[:Iqinj_lim][:min]
    I_qh1 = inner_controller_parameters[:Iqinj_lim][:max]
    V_ref0 = inner_controller_parameters[:V_ref0]
    T_iq = inner_controller_parameters[:T_iq]

    #Read local states
    Vt_filt = inner_controller_states[1]
    I_icv = inner_controller_states[2]

    #Compute additional states
    V_err = deadband_function(V_ref0 - Vt_filt, dbd1, dbd2)
    Iq_inj = clamp(K_qv * V_err, I_ql1, I_qh1)
    Iq_cmd = I_icv + Iq_inj
    Ip_min, Ip_max, Iq_min, Iq_max =
        current_limit_logic(inner_control, Val(PQ_Flag), Vt_filt, Ip_oc, Iq_cmd)
    Iq_cmd = clamp(Iq_cmd, Iq_min, Iq_max)
    Ip_cmd = clamp(Ip_oc, Ip_min, Ip_max)

    #ODE update
    inner_controller_ode[1] = low_pass_mass_matrix(V_t, Vt_filt, 1.0, T_rv)[2]
    inner_controller_ode[2] = low_pass_mass_matrix(Iq_oc, I_icv, 1.0, T_iq)[2]

    #Update Inner Vars
    inner_vars[Id_ic_var] = Ip_cmd
    inner_vars[Iq_ic_var] = Iq_cmd
    return
end

#Q_Flag = 1
function _mdl_ode_RE_inner_controller_B!(
    inner_controller_ode::AbstractArray{<:ACCEPTED_REAL_TYPES},
    inner_controller_states::AbstractArray{<:ACCEPTED_REAL_TYPES},
    inner_controller_parameters::AbstractArray{<:ACCEPTED_REAL_TYPES},
    ::Val{1},
    inner_control::PSY.RECurrentControlB,
    dynamic_device::DynamicWrapper{
        PSY.DynamicInverter{C, O, PSY.RECurrentControlB, DC, P, F, L},
    },
    inner_vars::AbstractVector,
    h,
    t,
) where {
    C <: PSY.Converter,
    O <: PSY.OuterControl,
    DC <: PSY.DCSource,
    P <: PSY.FrequencyEstimator,
    F <: PSY.Filter,
    L <: Union{Nothing, PSY.OutputCurrentLimiter},
}
    #Obtain inner variables for component
    V_t = sqrt(inner_vars[Vr_filter_var]^2 + inner_vars[Vi_filter_var]^2)
    Ip_oc = inner_vars[Id_oc_var]
    V_oc = inner_vars[V_oc_var]

    #Get Current Controller parameters
    PQ_Flag = PSY.get_PQ_Flag(inner_control)
    T_rv = inner_controller_parameters[:T_rv]
    dbd1 = inner_controller_parameters[:dbd_pnts1]
    dbd2 = inner_controller_parameters[:dbd_pnts2]
    K_qv = inner_controller_parameters[:K_qv]
    I_ql1 = inner_controller_parameters[:Iqinj_lim][:min]
    I_qh1 = inner_controller_parameters[:Iqinj_lim][:max]
    K_vp = inner_controller_parameters[:K_vp]
    K_vi = inner_controller_parameters[:K_vi]
    V_ref0 = inner_controller_parameters[:V_ref0]

    #Read local states
    Vt_filt = inner_controller_states[1]
    ξ_icv = inner_controller_states[2]

    #Compute additional states
    V_err = deadband_function(V_ref0 - Vt_filt, dbd1, dbd2)
    Iq_inj = clamp(K_qv * V_err, I_ql1, I_qh1)

    #Compute block derivatives
    _, dVtfilt_dt = low_pass_mass_matrix(V_t, Vt_filt, 1.0, T_rv)
    I_icv, dξicv_dt = pi_block(V_oc, ξ_icv, K_vp, K_vi)
    Iq_cmd = I_icv + Iq_inj
    Ip_min, Ip_max, Iq_min, Iq_max =
        current_limit_logic(inner_control, Val(PQ_Flag), Vt_filt, Ip_oc, Iq_cmd)
    Iq_cmd = clamp(Iq_cmd, Iq_min, Iq_max)
    Ip_cmd = clamp(Ip_oc, Ip_min, Ip_max)

    #ODE update
    inner_controller_ode[1] = dVtfilt_dt
    inner_controller_ode[2] = dξicv_dt

    #Update Inner Vars
    inner_vars[Id_ic_var] = Ip_cmd
    inner_vars[Iq_ic_var] = Iq_cmd
    return
end

############################################
### ODE calculations via device dispatch ###
############################################

function mdl_inner_ode!(
    device_states::AbstractArray{<:ACCEPTED_REAL_TYPES},
    output_ode::AbstractArray{<:ACCEPTED_REAL_TYPES},
    p::AbstractArray{<:ACCEPTED_REAL_TYPES},
    inner_vars::AbstractArray{<:ACCEPTED_REAL_TYPES},
    dynamic_device::DynamicWrapper{
        PSY.DynamicInverter{
            C,
            O,
            PSY.VoltageModeControl,
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
    O <: PSY.OuterControl,
    DC <: PSY.DCSource,
    P <: PSY.FrequencyEstimator,
    F <: PSY.Filter,
    L <: Union{Nothing, PSY.OutputCurrentLimiter},
}

    #Obtain external states inputs for component
    external_ix = get_input_port_ix(dynamic_device, PSY.VoltageModeControl)
    Ir_filter = device_states[external_ix[1]]
    Ii_filter = device_states[external_ix[2]]
    Ir_cnv = device_states[external_ix[3]]
    Ii_cnv = device_states[external_ix[4]]
    Vr_filter = device_states[external_ix[5]]
    Vi_filter = device_states[external_ix[6]]

    #Obtain inner variables for component
    ω_oc = inner_vars[ω_oc_var]
    θ_oc = inner_vars[θ_oc_var]
    v_refr = inner_vars[V_oc_var]
    Vdc = inner_vars[Vdc_var]

    cf = p[:params][:Filter][:cf]
    lf = p[:params][:Filter][:lf]

    #Get Voltage Controller parameters
    params = p[:params][:InnerControl]
    kpv = params[:kpv]
    kiv = params[:kiv]
    kffv = params[:kffv]
    rv = params[:rv]
    lv = params[:lv]
    kpc = params[:kpc]
    kic = params[:kic]
    kffi = params[:kffi]
    ωad = params[:ωad]
    kad = params[:kad]

    #Obtain indices for component w/r to device
    local_ix = get_local_state_ix(dynamic_device, PSY.VoltageModeControl)

    #Define internal states for frequency estimator
    internal_states = @view device_states[local_ix]
    ξ_d = internal_states[1]
    ξ_q = internal_states[2]
    γ_d = internal_states[3]
    γ_q = internal_states[4]
    ϕ_d = internal_states[5]
    ϕ_q = internal_states[6]

    #Transformations to dq frame
    I_dq_filter = ri_dq(θ_oc + pi / 2) * [Ir_filter; Ii_filter]
    I_dq_cnv = ri_dq(θ_oc + pi / 2) * [Ir_cnv; Ii_cnv]
    V_dq_filter = ri_dq(θ_oc + pi / 2) * [Vr_filter; Vi_filter]

    ### Compute 6 states ODEs (D'Arco EPSR122 Model) ###
    ## SRF Voltage Control w/ Virtual Impedance ##
    #Virtual Impedance
    Vd_filter_ref = (v_refr - rv * I_dq_filter[d] + ω_oc * lv * I_dq_filter[q])
    Vq_filter_ref = (-rv * I_dq_filter[q] - ω_oc * lv * I_dq_filter[d])

    #Voltage Control PI Blocks
    Id_pi, dξd_dt = pi_block(Vd_filter_ref - V_dq_filter[d], ξ_d, kpv, kiv)
    Iq_pi, dξq_dt = pi_block(Vq_filter_ref - V_dq_filter[q], ξ_q, kpv, kiv)
    #PI Integrator (internal state)
    output_ode[local_ix[1]] = dξd_dt
    output_ode[local_ix[2]] = dξq_dt

    #Compensate output Control Signal - Links to SRF Current Controller
    Id_cnv_ref = Id_pi - cf * ω_oc * V_dq_filter[q] + kffi * I_dq_filter[d]
    Iq_cnv_ref = Iq_pi + cf * ω_oc * V_dq_filter[d] + kffi * I_dq_filter[q]

    #Get limiter and apply output current limiting
    limiter = PSY.get_limiter(dynamic_device)
    Id_cnv_ref2, Iq_cnv_ref2 = limit_output_current(limiter, Id_cnv_ref, Iq_cnv_ref)

    ## SRF Current Control ##
    #Current Control PI Blocks
    Vd_pi, dγd_dt = pi_block(Id_cnv_ref2 - I_dq_cnv[d], γ_d, kpc, kic)
    Vq_pi, dγq_dt = pi_block(Iq_cnv_ref2 - I_dq_cnv[q], γ_q, kpc, kic)
    #PI Integrator (internal state)
    output_ode[local_ix[3]] = dγd_dt
    output_ode[local_ix[4]] = dγq_dt

    #Compensate References for Converter Output Voltage
    Vd_cnv_ref =
        Vd_pi - ω_oc * lf * I_dq_cnv[q] + kffv * V_dq_filter[d] -
        kad * (V_dq_filter[d] - ϕ_d)
    Vq_cnv_ref =
        Vq_pi + ω_oc * lf * I_dq_cnv[d] + kffv * V_dq_filter[q] -
        kad * (V_dq_filter[q] - ϕ_q)

    #Active Damping LPF (internal state)
    output_ode[local_ix[5]] = low_pass(V_dq_filter[d], ϕ_d, 1.0, 1.0 / ωad)[2]
    output_ode[local_ix[6]] = low_pass(V_dq_filter[q], ϕ_q, 1.0, 1.0 / ωad)[2]

    #Update inner_vars
    #Modulation Commands to Converter
    inner_vars[md_var] = Vd_cnv_ref / Vdc
    inner_vars[mq_var] = Vq_cnv_ref / Vdc
    return
end

function mdl_inner_ode!(
    device_states::AbstractArray{<:ACCEPTED_REAL_TYPES},
    output_ode::AbstractArray{<:ACCEPTED_REAL_TYPES},
    inner_vars::AbstractArray{<:ACCEPTED_REAL_TYPES},
    dynamic_device::DynamicWrapper{
        PSY.DynamicInverter{
            C,
            O,
            PSY.VoltageModeControl,
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
    O <: PSY.OuterControl,
    DC <: PSY.DCSource,
    P <: PSY.FrequencyEstimator,
    F <: PSY.Filter,
    L <: Union{PSY.HybridOutputCurrentLimiter, PSY.SaturationOutputCurrentLimiter},
}

    #Obtain external states inputs for component
    external_ix = get_input_port_ix(dynamic_device, PSY.VoltageModeControl)
    Ir_filter = device_states[external_ix[1]]
    Ii_filter = device_states[external_ix[2]]
    Ir_cnv = device_states[external_ix[3]]
    Ii_cnv = device_states[external_ix[4]]
    Vr_filter = device_states[external_ix[5]]
    Vi_filter = device_states[external_ix[6]]

    #Obtain inner variables for component
    ω_oc = inner_vars[ω_oc_var]
    θ_oc = inner_vars[θ_oc_var]
    v_refr = inner_vars[V_oc_var]
    Vdc = inner_vars[Vdc_var]

    #Get Voltage Controller parameters
    inner_control = PSY.get_inner_control(dynamic_device)
    filter = PSY.get_filter(dynamic_device)
    kpv = PSY.get_kpv(inner_control)
    kiv = PSY.get_kiv(inner_control)
    kffi = PSY.get_kffi(inner_control)
    cf = PSY.get_cf(filter)
    rv = PSY.get_rv(inner_control)
    lv = PSY.get_lv(inner_control)

    #Get Current Controller parameters
    kpc = PSY.get_kpc(inner_control)
    kic = PSY.get_kic(inner_control)
    kffv = PSY.get_kffv(inner_control)
    lf = PSY.get_lf(filter)
    ωad = PSY.get_ωad(inner_control)
    kad = PSY.get_kad(inner_control)

    #Obtain indices for component w/r to device
    local_ix = get_local_state_ix(dynamic_device, PSY.VoltageModeControl)

    #Define internal states for frequency estimator
    internal_states = @view device_states[local_ix]
    ξ_d = internal_states[1]
    ξ_q = internal_states[2]
    γ_d = internal_states[3]
    γ_q = internal_states[4]
    ϕ_d = internal_states[5]
    ϕ_q = internal_states[6]

    #Transformations to dq frame
    I_dq_filter = ri_dq(θ_oc + pi / 2) * [Ir_filter; Ii_filter]
    I_dq_cnv = ri_dq(θ_oc + pi / 2) * [Ir_cnv; Ii_cnv]
    V_dq_filter = ri_dq(θ_oc + pi / 2) * [Vr_filter; Vi_filter]

    ### Compute 6 states ODEs (D'Arco EPSR122 Model) ###
    ## SRF Voltage Control w/ Virtual Impedance ##
    #Virtual Impedance
    Vd_filter_ref = (v_refr - rv * I_dq_filter[d] + ω_oc * lv * I_dq_filter[q])
    Vq_filter_ref = (-rv * I_dq_filter[q] - ω_oc * lv * I_dq_filter[d])

    #Voltage Control PI Blocks
    Prop_d = kpv * (Vd_filter_ref - V_dq_filter[d])
    Prop_q = kpv * (Vq_filter_ref - V_dq_filter[q])
    Integral_d = ξ_d
    Integral_q = ξ_q

    #Compensate output Control Signal - Links to SRF Current Controller
    Id_cnv_ref =
        Prop_d + (kiv * Integral_d) - cf * ω_oc * V_dq_filter[q] + kffi * I_dq_filter[d]
    Iq_cnv_ref =
        Prop_q + (kiv * Integral_q) + cf * ω_oc * V_dq_filter[d] + kffi * I_dq_filter[q]

    #Get limiter and apply output current limiting
    limiter = PSY.get_limiter(dynamic_device)
    Id_cnv_ref2, Iq_cnv_ref2, Del_Vv_d, Del_Vv_q =
        limit_output_current(limiter, Id_cnv_ref, Iq_cnv_ref, ω_oc)

    # Limiter anti- windup 
    dξd_dt = (Vd_filter_ref - V_dq_filter[d]) - Del_Vv_d
    dξq_dt = (Vq_filter_ref - V_dq_filter[q]) - Del_Vv_q

    #PI Integrator (internal state) 
    output_ode[local_ix[1]] = dξd_dt
    output_ode[local_ix[2]] = dξq_dt

    ## SRF Current Control ##
    #Current Control PI Blocks
    Vd_pi, dγd_dt = pi_block(Id_cnv_ref2 - I_dq_cnv[d], γ_d, kpc, kic)
    Vq_pi, dγq_dt = pi_block(Iq_cnv_ref2 - I_dq_cnv[q], γ_q, kpc, kic)
    #PI Integrator (internal state)
    output_ode[local_ix[3]] = dγd_dt
    output_ode[local_ix[4]] = dγq_dt

    #Compensate References for Converter Output Voltage
    Vd_cnv_ref =
        Vd_pi - ω_oc * lf * I_dq_cnv[q] + kffv * V_dq_filter[d] -
        kad * (V_dq_filter[d] - ϕ_d)
    Vq_cnv_ref =
        Vq_pi + ω_oc * lf * I_dq_cnv[d] + kffv * V_dq_filter[q] -
        kad * (V_dq_filter[q] - ϕ_q)

    #Active Damping LPF (internal state)
    output_ode[local_ix[5]] = low_pass(V_dq_filter[d], ϕ_d, 1.0, 1.0 / ωad)[2]
    output_ode[local_ix[6]] = low_pass(V_dq_filter[q], ϕ_q, 1.0, 1.0 / ωad)[2]

    #Update inner_vars
    #Modulation Commands to Converter
    inner_vars[md_var] = Vd_cnv_ref / Vdc
    inner_vars[mq_var] = Vq_cnv_ref / Vdc
    return
end

function mdl_inner_ode!(
    device_states::AbstractArray{<:ACCEPTED_REAL_TYPES},
    output_ode::AbstractArray{<:ACCEPTED_REAL_TYPES},
    p::AbstractArray{<:ACCEPTED_REAL_TYPES},
    inner_vars::AbstractArray{<:ACCEPTED_REAL_TYPES},
    dynamic_device::DynamicWrapper{
        PSY.DynamicInverter{C, O, PSY.CurrentModeControl, DC, P, F, L},
    },
    h,
    t,
) where {
    C <: PSY.Converter,
    O <: PSY.OuterControl,
    DC <: PSY.DCSource,
    P <: PSY.FrequencyEstimator,
    F <: PSY.Filter,
    L <: Union{Nothing, PSY.OutputCurrentLimiter},
}

    #Obtain external states inputs for component
    external_ix = get_input_port_ix(dynamic_device, PSY.CurrentModeControl)
    # Ir_filter = device_states[external_ix[1]]
    # Ii_filter = device_states[external_ix[2]]
    Ir_cnv = device_states[external_ix[3]]
    Ii_cnv = device_states[external_ix[4]]
    Vr_filter = device_states[external_ix[5]]
    Vi_filter = device_states[external_ix[6]]

    #Obtain inner variables for component
    ω_oc = inner_vars[ω_oc_var]
    θ_oc = inner_vars[θ_oc_var]
    Vdc = inner_vars[Vdc_var]

    #Get Current Controller parameters
    kpc = p[:params][:InnerControl][:kpc]
    kic = p[:params][:InnerControl][:kic]
    kffv = p[:params][:InnerControl][:kffv]
    lf = p[:params][:Filter][:lf]

    #Obtain indices for component w/r to device
    local_ix = get_local_state_ix(dynamic_device, PSY.CurrentModeControl)

    #Define internal states for frequency estimator
    internal_states = @view device_states[local_ix]
    γ_d = internal_states[1]
    γ_q = internal_states[2]

    #Transformations to dq frame
    I_dq_cnv = ri_dq(θ_oc + pi / 2) * [Ir_cnv; Ii_cnv]
    V_dq_filter = ri_dq(θ_oc + pi / 2) * [Vr_filter; Vi_filter]

    #Input Control Signal - Links to SRF Current Controller
    Id_cnv_ref = inner_vars[Id_oc_var]
    Iq_cnv_ref = inner_vars[Iq_oc_var]

    #Current Control PI Blocks
    Vd_pi, dγd_dt = pi_block(Id_cnv_ref - I_dq_cnv[d], γ_d, kpc, kic)
    Vq_pi, dγq_dt = pi_block(Iq_cnv_ref - I_dq_cnv[q], γ_q, kpc, kic)
    #PI Integrator (internal state)
    output_ode[local_ix[1]] = dγd_dt
    output_ode[local_ix[2]] = dγq_dt

    #Compensate references for Converter Output Voltage
    Vd_cnv_ref = Vd_pi - ω_oc * lf * I_dq_cnv[q] + kffv * V_dq_filter[d]
    Vq_cnv_ref = Vq_pi + ω_oc * lf * I_dq_cnv[d] + kffv * V_dq_filter[q]

    #Update inner_vars
    #Modulation Commands to Converter
    inner_vars[md_var] = Vd_cnv_ref / Vdc
    inner_vars[mq_var] = Vq_cnv_ref / Vdc
    return
end

function mdl_inner_ode!(
    device_states::AbstractArray{<:ACCEPTED_REAL_TYPES},
    output_ode::AbstractArray{<:ACCEPTED_REAL_TYPES},
    p::AbstractArray{<:ACCEPTED_REAL_TYPES},
    inner_vars::AbstractArray{<:ACCEPTED_REAL_TYPES},
    dynamic_device::DynamicWrapper{
        PSY.DynamicInverter{C, O, PSY.RECurrentControlB, DC, P, F, L},
    },
    h,
    t,
) where {
    C <: PSY.Converter,
    O <: PSY.OuterControl,
    DC <: PSY.DCSource,
    P <: PSY.FrequencyEstimator,
    F <: PSY.Filter,
    L <: Union{Nothing, PSY.OutputCurrentLimiter},
}
    #Get Current Controller parameters
    inner_control = PSY.get_inner_control(dynamic_device)
    Q_Flag = PSY.get_Q_Flag(inner_control)

    #Obtain indices for component w/r to device
    local_ix = get_local_state_ix(dynamic_device, PSY.RECurrentControlB)
    #Define internal states for Inner Control
    internal_states = @view device_states[local_ix]
    internal_ode = @view output_ode[local_ix]
    internal_parameters = @view p[:params][:InnerControl]
    # TODO: Voltage Dip Freeze logic

    #Dispatch inner controller ODE calculation
    _mdl_ode_RE_inner_controller_B!(
        internal_ode,
        internal_states,
        internal_parameters,
        Val(Q_Flag),
        inner_control,
        dynamic_device,
        inner_vars,
        h,
        t,
    )
    return
end
