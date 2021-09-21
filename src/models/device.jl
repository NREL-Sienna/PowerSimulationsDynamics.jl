function device_mass_matrix_entries!(::AbstractArray, ::DynamicWrapper{T}) where {T}
    error("Mass Matrix not implemented for models $T")
end

function device_mass_matrix_entries!(
    mass_matrix::AbstractArray,
    dynamic_device::DynamicWrapper{DynG},
) where {DynG <: PSY.DynamicGenerator}
    global_index = get_global_index(dynamic_device)
    mass_matrix_tg_entries!(mass_matrix, PSY.get_prime_mover(dynamic_device), global_index)
    mass_matrix_pss_entries!(mass_matrix, PSY.get_pss(dynamic_device), global_index)
    mass_matrix_avr_entries!(mass_matrix, PSY.get_avr(dynamic_device), global_index)
    mass_matrix_machine_entries!(mass_matrix, PSY.get_machine(dynamic_device), global_index)
    mass_matrix_shaft_entries!(mass_matrix, PSY.get_shaft(dynamic_device), global_index)
    return
end

function device!(
    device_states::AbstractArray{T},
    output_ode::AbstractArray{T},
    voltage_r::T,
    voltage_i::T,
    current_r::AbstractArray{T},
    current_i::AbstractArray{T},
    global_vars::AbstractArray{T},
    inner_vars::AbstractArray{T},
    dynamic_device::DynamicWrapper{DynG},
    t,
) where {DynG <: PSY.DynamicGenerator, T <: Real}
    inner_vars[VR_gen_var] = voltage_r
    inner_vars[VI_gen_var] = voltage_i

    sys_ω = global_vars[GLOBAL_VAR_SYS_FREQ_INDEX]

    #Obtain ODEs and Mechanical Power for Turbine Governor
    mdl_tg_ode!(device_states, output_ode, inner_vars, sys_ω, dynamic_device)

    #Obtain ODEs for AVR
    mdl_pss_ode!(device_states, output_ode, inner_vars, sys_ω, dynamic_device)

    #Obtain ODEs for AVR
    mdl_avr_ode!(device_states, output_ode, inner_vars, dynamic_device)

    #Obtain ODEs for Machine
    mdl_machine_ode!(
        device_states,
        output_ode,
        inner_vars,
        current_r,
        current_i,
        dynamic_device,
    )

    #Obtain ODEs for PSY.Shaft
    mdl_shaft_ode!(device_states, output_ode, inner_vars, sys_ω, dynamic_device)
    return
end

function device!(
    voltage_r::T,
    voltage_i::T,
    current_r::AbstractArray{T},
    current_i::AbstractArray{T},
    ::AbstractArray{T},
    ::AbstractArray{T},
    device::StaticWrapper{PSY.Source, U},
    t,
) where {T <: Real, U <: BusCategory}
    mdl_source!(voltage_r, voltage_i, current_r, current_i, device)
    return
end

function device!(
    voltage_r::T,
    voltage_i::T,
    current_r::AbstractArray{T},
    current_i::AbstractArray{T},
    ::AbstractArray{T},
    ::AbstractArray{T},
    device::StaticWrapper{PSY.PowerLoad, U},
    t,
) where {T <: Real, U <: LoadCategory}
    # For now model all loads as constant impedance
    mdl_Zload!(voltage_r, voltage_i, current_r, current_i, device.device)
    return
end

function device_mass_matrix_entries!(
    mass_matrix::AbstractArray{Float64},
    dynamic_device::DynamicWrapper{DynI},
) where {DynI <: PSY.DynamicInverter}
    global_index = get_global_index(dynamic_device)
    mass_matrix_DCside_entries!(
        mass_matrix,
        PSY.get_dc_source(dynamic_device),
        global_index,
    )
    mass_matrix_freq_estimator_entries!(
        mass_matrix,
        PSY.get_freq_estimator(dynamic_device),
        global_index,
    )
    mass_matrix_outer_entries!(
        mass_matrix,
        PSY.get_outer_control(dynamic_device),
        global_index,
    )
    mass_matrix_inner_entries!(
        mass_matrix,
        PSY.get_inner_control(dynamic_device),
        global_index,
    )
    mass_matrix_converter_entries!(
        mass_matrix,
        PSY.get_converter(dynamic_device),
        global_index,
    )
    mass_matrix_filter_entries!(mass_matrix, PSY.get_filter(dynamic_device), global_index)
    return
end

function device!(
    device_states::AbstractArray{T},
    output_ode::AbstractArray{T},
    voltage_r::T,
    voltage_i::T,
    current_r::AbstractArray{T},
    current_i::AbstractArray{T},
    global_vars::AbstractArray{T},
    inner_vars::AbstractArray{T},
    dynamic_device::DynamicWrapper{DynI},
    t,
) where {DynI <: PSY.DynamicInverter, T <: Real}
    #Obtain global vars
    sys_ω = global_vars[GLOBAL_VAR_SYS_FREQ_INDEX]

    #Update Voltage data
    inner_vars[Vr_inv_var] = voltage_r
    inner_vars[Vi_inv_var] = voltage_i

    #Update V_ref
    V_ref = get_V_ref(dynamic_device)
    inner_vars[V_oc_var] = V_ref

    #Update current inner_vars
    _update_inner_vars!(device_states, output_ode, sys_ω, inner_vars, dynamic_device)

    #Obtain ODES for DC side
    mdl_DCside_ode!(device_states, output_ode, sys_ω, inner_vars, dynamic_device)

    #Obtain ODEs for PLL
    mdl_freq_estimator_ode!(device_states, output_ode, inner_vars, sys_ω, dynamic_device)

    #Obtain ODEs for OuterLoop
    mdl_outer_ode!(device_states, output_ode, inner_vars, sys_ω, dynamic_device)

    #Obtain inner controller ODEs and modulation commands
    mdl_inner_ode!(device_states, output_ode, inner_vars, dynamic_device)

    #Obtain converter relations
    mdl_converter_ode!(device_states, output_ode, inner_vars, dynamic_device)

    #Obtain ODEs for output filter
    mdl_filter_ode!(
        device_states,
        output_ode,
        current_r,
        current_i,
        inner_vars,
        sys_ω,
        dynamic_device,
    )

    return
end

function device_mass_matrix_entries!(
    mass_matrix::AbstractArray,
    dynamic_device::DynamicWrapper{PSY.PeriodicVariableSource},
)
    global_index = get_global_index(dynamic_device)
    mass_matrix_pvs_entries!(mass_matrix, dynamic_device, global_index)
    return
end

function mass_matrix_pvs_entries!(
    mass_matrix,
    pvs::DynamicWrapper{PSY.PeriodicVariableSource},
    global_index::ImmutableDict{Symbol, Int64},
)
    @debug "Using default mass matrix entries $pvs"
end

function device!(
    device_states::AbstractArray{T},
    output_ode::AbstractArray{T},
    voltage_r::AbstractArray{T},
    voltage_i::AbstractArray{T},
    current_r::AbstractArray{T},
    current_i::AbstractArray{T},
    ::AbstractArray{T},
    dynamic_device::PSY.PeriodicVariableSource,
    ::AbstractArray{T},
    t,
) where {T <: Real}
    ω_θ = PSY.get_internal_angle_frequencies(dynamic_device)
    ω_V = PSY.get_internal_angle_frequencies(dynamic_device)

    dV = 0
    for (ix, A) in enumerate(PSY.get_internal_voltage_coefficients(dynamic_device))
        t <= 0 && continue
        dV += ω_V[ix] * (A[1] * cos(ω_V[ix] * t) - A[2] * sin(ω_V[ix] * t))
    end

    dθ = 0
    for (ix, A) in enumerate(PSY.get_internal_angle_coefficients(dynamic_device))
        t <= 0 && continue
        dθ += ω_θ[ix] * (A[1] * cos(ω_θ[ix] * t) - A[2] * sin(ω_θ[ix] * t))
    end

    # Internal Voltage states
    V_R = device_states[1] * cos(device_states[2])
    V_I = device_states[1] * sin(device_states[2])
    output_ode[1] = dV
    output_ode[2] = dθ

    #update current
    R_th = PSY.get_R_th(dynamic_device)
    X_th = PSY.get_X_th(dynamic_device)
    Zmag = R_th^2 + X_th^2
    current_r[1] += R_th * (V_R - voltage_r[1]) / Zmag + X_th * (V_I - voltage_i[1]) / Zmag #in system pu flowing out
    current_i[1] += R_th * (V_I - voltage_i[1]) / Zmag - X_th * (V_R - voltage_r[1]) / Zmag #in system pu flowing out

    return
end

function _update_inner_vars!(
    device_states,
    output_ode,
    ω_sys,
    inner_vars,
    dynamic_device::DynamicWrapper{PSY.DynamicInverter{C, O, IC, DC, P, F}},
) where {
    C <: PSY.Converter,
    O <: PSY.OuterControl,
    IC <: PSY.InnerControl,
    DC <: PSY.DCSource,
    P <: PSY.FrequencyEstimator,
    F <: PSY.Filter,
}
    #do nothing
end

function _update_inner_vars!(
    device_states,
    output_ode,
    ω_sys,
    inner_vars,
    dynamic_device::DynamicWrapper{
        PSY.DynamicInverter{PSY.RenewableEnergyConverterTypeA, O, IC, DC, P, PSY.RLFilter},
    },
) where {
    O <: PSY.OuterControl,
    IC <: PSY.InnerControl,
    DC <: PSY.DCSource,
    P <: PSY.FrequencyEstimator,
}
    V_R = inner_vars[Vr_inv_var]
    V_I = inner_vars[Vi_inv_var]
    V_t = sqrt(V_R^2 + V_I^2)
    θ = atan(V_I / V_R)

    #Get Converter parameters
    converter = PSY.get_converter(dynamic_device)
    Brkpt = PSY.get_Brkpt(converter)
    Zerox = PSY.get_Zerox(converter)
    Lvpl1 = PSY.get_Lvpl1(converter)
    Vo_lim = PSY.get_Vo_lim(converter)
    Lv_pnt0, Lv_pnt1 = PSY.get_Lv_pnts(converter)
    K_hv = PSY.get_K_hv(converter)
    Lvpl_sw = PSY.get_Lvpl_sw(converter)
    R_source = PSY.get_R_source(converter)
    X_source = PSY.get_X_source(converter)

    #Define internal states for Converter
    converter_ix = get_local_state_ix(dynamic_device, PSY.RenewableEnergyConverterTypeA)
    converter_states = @view device_states[converter_ix]
    Ip = converter_states[1]
    Iq = converter_states[2]
    Vmeas = converter_states[3]

    #Saturate Ip if LVPL is active
    Ip_sat = Ip
    if Lvpl_sw == 1
        LVPL = get_LVPL_gain(Vmeas, Zerox, Brkpt, Lvpl1)
        Ip_sat = Ip <= LVPL ? Ip : LVPL
    end
    G_lv = get_LV_current_gain(V_t, Lv_pnt0, Lv_pnt1)
    Iq_extra = max(K_hv * (V_t - Vo_lim), 0.0)
    Id_cnv = G_lv * Ip_sat
    Iq_cnv = -Iq - Iq_extra
    #Reference Transformation
    Ir_cnv = Id_cnv * cos(θ) - Iq_cnv * sin(θ)
    Ii_cnv = Id_cnv * sin(θ) + Iq_cnv * cos(θ)

    #Obtain parameters
    filt = PSY.get_filter(dynamic_device)
    rf = PSY.get_rf(filt)
    lf = PSY.get_lf(filt)

    function V_cnv_calc(Ir_cnv, Ii_cnv, Vr_inv, Vi_inv)
        if lf != 0.0 || rf != 0.0
            Z_source_mag_sq = R_source^2 + X_source^2
            Zf_mag_sq = rf^2 + lf^2
            r_total_ratio = rf / Zf_mag_sq + R_source / Z_source_mag_sq
            l_total_ratio = lf / Zf_mag_sq + X_source / Z_source_mag_sq
            denom = 1.0 / (r_total_ratio^2 + l_total_ratio^2)
            rf_ratio = rf / Zf_mag_sq
            denom * (
                (Ir_cnv + Vi_inv * lf_ratio + Vr_inv * rf_ratio) * r_total_ratio -
                (Ii_cnv + Vi_inv * rf_ratio - Vr_inv * lf_ratio) * l_total_ratio
            )
            Vi_cnv =
                denom * (
                    (Ir_cnv + Vi_inv * lf_ratio + Vr_inv * rf_ratio) * l_total_ratio +
                    (Ii_cnv + Vi_inv * rf_ratio - Vr_inv * lf_ratio) * r_total_ratio
                )
        else
            Vr_cnv = Vr_inv
            Vi_cnv = Vi_inv
        end
        return Vr_cnv, Vi_cnv
    end

    Vr_cnv, Vi_cnv = V_cnv_calc(Ir_cnv, Ii_cnv, V_R, V_I)

    #Compute output currents
    if lf != 0.0 || rf != 0.0
        Vr_inv = V_R
        Vi_inv = V_I
        Zmag_squared = rf^2 + lf^2
        Ir_filt = (1.0 / Zmag_squared) * ((Vr_cnv - Vr_inv) * rf + (Vi_cnv - Vi_inv) * lf)
        Ii_filt = (1.0 / Zmag_squared) * ((Vi_cnv - Vi_inv) * rf - (Vr_cnv - Vr_inv) * lf)
    else
        Ir_filt = Ir_cnv
        Ii_filt = Ii_cnv
    end

    #Update inner_vars
    inner_vars[Ir_cnv_var] = Ir_cnv
    inner_vars[Ii_cnv_var] = Ii_cnv
    inner_vars[Vr_cnv_var] = Vr_cnv
    inner_vars[Vi_cnv_var] = Vi_cnv
    inner_vars[Ir_inv_var] = Ir_filt
    inner_vars[Ii_inv_var] = Ii_filt
end
