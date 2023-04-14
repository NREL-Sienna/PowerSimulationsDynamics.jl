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
) where {DynG <: PSY.DynamicGenerator, T <: ACCEPTED_REAL_TYPES}
    if get_connection_status(dynamic_device) < 1.0
        output_ode .= zero(T)
        return
    end

    inner_vars[VR_gen_var] = voltage_r
    inner_vars[VI_gen_var] = voltage_i

    sys_ω = global_vars[GLOBAL_VAR_SYS_FREQ_INDEX]

    #Update Inner Vars
    _update_inner_vars!(device_states, output_ode, sys_ω, inner_vars, dynamic_device)

    #Obtain ODEs and Mechanical Power for Turbine Governor
    mdl_tg_ode!(device_states, output_ode, inner_vars, sys_ω, dynamic_device)

    #Obtain ODEs for PSS
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
) where {T <: ACCEPTED_REAL_TYPES, U <: BusCategory}
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
    device::StaticLoadWrapper,
    t,
) where {T <: ACCEPTED_REAL_TYPES}
    mdl_zip_load!(voltage_r, voltage_i, current_r, current_i, device)
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
    f0 = get_system_base_frequency(dynamic_device)
    mass_matrix_filter_entries!(
        mass_matrix,
        PSY.get_filter(dynamic_device),
        global_index,
        f0,
    )
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
) where {DynI <: PSY.DynamicInverter, T <: ACCEPTED_REAL_TYPES}
    if get_connection_status(dynamic_device) < 1.0
        output_ode .= zero(T)
        return
    end

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
    voltage_r::T,
    voltage_i::T,
    current_r::AbstractArray{T},
    current_i::AbstractArray{T},
    ::AbstractArray{T},
    ::AbstractArray{T},
    dynamic_device::DynamicWrapper{PSY.PeriodicVariableSource},
    t,
) where {T <: ACCEPTED_REAL_TYPES}
    ω_θ = PSY.get_internal_angle_frequencies(get_device(dynamic_device))
    ω_V = PSY.get_internal_angle_frequencies(get_device(dynamic_device))

    dV = 0
    for (ix, A) in
        enumerate(PSY.get_internal_voltage_coefficients(get_device(dynamic_device)))
        t <= 0 && continue
        dV += ω_V[ix] * (A[1] * cos(ω_V[ix] * t) - A[2] * sin(ω_V[ix] * t))
    end

    dθ = 0
    for (ix, A) in
        enumerate(PSY.get_internal_angle_coefficients(get_device(dynamic_device)))
        t <= 0 && continue
        dθ += ω_θ[ix] * (A[1] * cos(ω_θ[ix] * t) - A[2] * sin(ω_θ[ix] * t))
    end

    # Internal Voltage states
    V_R = device_states[1] * cos(device_states[2])
    V_I = device_states[1] * sin(device_states[2])
    output_ode[1] = dV
    output_ode[2] = dθ

    #update current
    R_th = PSY.get_R_th(get_device(dynamic_device))
    X_th = PSY.get_X_th(get_device(dynamic_device))
    Zmag = R_th^2 + X_th^2
    current_r[1] += R_th * (V_R - voltage_r[1]) / Zmag + X_th * (V_I - voltage_i[1]) / Zmag #in system pu flowing out
    current_i[1] += R_th * (V_I - voltage_i[1]) / Zmag - X_th * (V_R - voltage_r[1]) / Zmag #in system pu flowing out

    return
end

function _update_inner_vars!(
    ::AbstractArray{<:ACCEPTED_REAL_TYPES},
    ::AbstractArray{<:ACCEPTED_REAL_TYPES},
    ω_sys::ACCEPTED_REAL_TYPES,
    ::AbstractArray{<:ACCEPTED_REAL_TYPES},
    dynamic_device::DynamicWrapper{PSY.DynamicGenerator{M, S, A, TG, P}},
) where {M <: PSY.Machine, S <: PSY.Shaft, A <: PSY.AVR, TG <: PSY.TurbineGov, P <: PSY.PSS}
    return
end

function _update_inner_vars!(
    device_states::AbstractArray{<:ACCEPTED_REAL_TYPES},
    ::AbstractArray{<:ACCEPTED_REAL_TYPES},
    ω_sys::ACCEPTED_REAL_TYPES,
    inner_vars::AbstractArray{<:ACCEPTED_REAL_TYPES},
    dynamic_device::DynamicWrapper{PSY.DynamicGenerator{M, S, A, TG, P}},
) where {
    M <: Union{PSY.RoundRotorQuadratic, PSY.RoundRotorExponential},
    S <: PSY.Shaft,
    A <: PSY.AVR,
    TG <: PSY.TurbineGov,
    P <: PSY.PSS,
}
    #Obtain indices for component w/r to device
    machine = PSY.get_machine(dynamic_device)
    local_ix = get_local_state_ix(dynamic_device, M)

    #Define internal states for component
    internal_states = @view device_states[local_ix]
    eq_p = internal_states[1]
    ed_p = internal_states[2]
    ψ_kd = internal_states[3]
    ψ_kq = internal_states[4]

    #Obtain external states inputs for component
    external_ix = get_input_port_ix(dynamic_device, M)
    δ = device_states[external_ix[1]]

    #Obtain inner variables for component
    V_tR = inner_vars[VR_gen_var]
    V_tI = inner_vars[VI_gen_var]

    #Get parameters
    R = PSY.get_R(machine)
    Xd = PSY.get_Xd(machine)
    Xd_p = PSY.get_Xd_p(machine)
    Xd_pp = PSY.get_Xd_pp(machine)
    Xq_pp = Xd_pp
    Xl = PSY.get_Xl(machine)
    γ_d1 = PSY.get_γ_d1(machine)
    γ_q1 = PSY.get_γ_q1(machine)
    γ_d2 = PSY.get_γ_d2(machine)

    #RI to dq transformation
    V_dq = ri_dq(δ) * [V_tR; V_tI]

    #Additional Fluxes
    ψq_pp = γ_q1 * ed_p + ψ_kq * (1 - γ_q1)
    ψd_pp = γ_d1 * eq_p + γ_d2 * (Xd_p - Xl) * ψ_kd
    ψ_pp = sqrt(ψd_pp^2 + ψq_pp^2)
    #Currents
    I_d =
        (1.0 / (R^2 + Xq_pp * Xd_pp)) *
        (-R * (V_dq[1] - ψq_pp) + Xq_pp * (-V_dq[2] + ψd_pp))
    Se = saturation_function(machine, ψ_pp)
    Xad_Ifd = eq_p + (Xd - Xd_p) * (γ_d1 * I_d - γ_d2 * ψ_kd + γ_d2 * eq_p) + Se * ψd_pp

    #Update Xad_Ifd
    inner_vars[Xad_Ifd_var] = Xad_Ifd
    return
end

function _update_inner_vars!(
    device_states::AbstractArray{<:ACCEPTED_REAL_TYPES},
    ::AbstractArray{<:ACCEPTED_REAL_TYPES},
    ω_sys::ACCEPTED_REAL_TYPES,
    inner_vars::AbstractArray{<:ACCEPTED_REAL_TYPES},
    dynamic_device::DynamicWrapper{
        PSY.DynamicGenerator{PSY.SalientPoleQuadratic, S, A, TG, P},
    },
) where {S <: PSY.Shaft, A <: PSY.AVR, TG <: PSY.TurbineGov, P <: PSY.PSS}
    #Obtain indices for component w/r to device
    machine = PSY.get_machine(dynamic_device)
    local_ix = get_local_state_ix(dynamic_device, typeof(machine))

    #Define internal states for component
    internal_states = @view device_states[local_ix]
    eq_p = internal_states[1]
    ψ_kd = internal_states[2]
    ψq_pp = internal_states[3]

    #Obtain external states inputs for component
    external_ix = get_input_port_ix(dynamic_device, typeof(machine))
    δ = device_states[external_ix[1]]

    #Obtain inner variables for component
    V_tR = inner_vars[VR_gen_var]
    V_tI = inner_vars[VI_gen_var]

    #Get parameters
    R = PSY.get_R(machine)
    Xd = PSY.get_Xd(machine)
    Xd_p = PSY.get_Xd_p(machine)
    Xd_pp = PSY.get_Xd_pp(machine)
    Xl = PSY.get_Xl(machine)
    γ_d1 = PSY.get_γ_d1(machine)
    γ_q1 = PSY.get_γ_q1(machine)
    γ_d2 = PSY.get_γ_d2(machine)

    #RI to dq transformation
    V_d, V_q = ri_dq(δ) * [V_tR; V_tI]

    #Additional Fluxes
    ψd_pp = γ_d1 * eq_p + γ_q1 * ψ_kd

    #Currents
    I_d = (1.0 / (R^2 + Xd_pp^2)) * (-R * (V_d + ψq_pp) + Xd_pp * (ψd_pp - V_q))
    Se = saturation_function(machine, eq_p)
    Xad_Ifd =
        eq_p + Se * eq_p + (Xd - Xd_p) * (I_d + γ_d2 * (eq_p - ψ_kd - (Xd_p - Xl) * I_d))

    #Update Xad_Ifd
    inner_vars[Xad_Ifd_var] = Xad_Ifd
    return
end

function _update_inner_vars!(
    device_states::AbstractArray{<:ACCEPTED_REAL_TYPES},
    ::AbstractArray{<:ACCEPTED_REAL_TYPES},
    ω_sys::ACCEPTED_REAL_TYPES,
    inner_vars::AbstractArray{<:ACCEPTED_REAL_TYPES},
    dynamic_device::DynamicWrapper{
        PSY.DynamicGenerator{PSY.SalientPoleExponential, S, A, TG, P},
    },
) where {S <: PSY.Shaft, A <: PSY.AVR, TG <: PSY.TurbineGov, P <: PSY.PSS}
    #Obtain indices for component w/r to device
    machine = PSY.get_machine(dynamic_device)
    local_ix = get_local_state_ix(dynamic_device, typeof(machine))

    #Define internal states for component
    internal_states = @view device_states[local_ix]
    eq_p = internal_states[1]
    ψ_kd = internal_states[2]
    ψq_pp = internal_states[3]

    #Obtain external states inputs for component
    external_ix = get_input_port_ix(dynamic_device, typeof(machine))
    δ = device_states[external_ix[1]]

    #Obtain inner variables for component
    V_tR = inner_vars[VR_gen_var]
    V_tI = inner_vars[VI_gen_var]

    #Get parameters
    R = PSY.get_R(machine)
    Xd = PSY.get_Xd(machine)
    Xd_p = PSY.get_Xd_p(machine)
    Xd_pp = PSY.get_Xd_pp(machine)
    Xq_pp = Xd_pp
    Xl = PSY.get_Xl(machine)
    γ_d1 = PSY.get_γ_d1(machine)
    γ_q1 = PSY.get_γ_q1(machine)
    γ_d2 = PSY.get_γ_d2(machine)

    #RI to dq transformation
    V_d, V_q = ri_dq(δ) * [V_tR; V_tI]

    #Additional Fluxes
    ψd_pp = γ_d1 * eq_p + γ_q1 * ψ_kd
    ψ_pp = sqrt(ψd_pp^2 + ψq_pp^2)

    #Currents
    I_d = (1.0 / (R^2 + Xd_pp^2)) * (-R * (V_d - ψq_pp) + Xq_pp * (-V_q + ψd_pp))
    Se = saturation_function(machine, ψ_pp)
    Xad_Ifd =
        eq_p + Se * ψd_pp + (Xd - Xd_p) * (I_d + γ_d2 * (eq_p - ψ_kd - (Xd_p - Xl) * I_d))

    #Update Xad_Ifd
    inner_vars[Xad_Ifd_var] = Xad_Ifd
    return
end

function _update_inner_vars!(
    ::AbstractArray{<:ACCEPTED_REAL_TYPES},
    ::AbstractArray{<:ACCEPTED_REAL_TYPES},
    ω_sys::ACCEPTED_REAL_TYPES,
    ::AbstractArray{<:ACCEPTED_REAL_TYPES},
    dynamic_device::DynamicWrapper{PSY.DynamicInverter{C, O, IC, DC, P, F}},
) where {
    C <: PSY.Converter,
    O <: PSY.OuterControl,
    IC <: PSY.InnerControl,
    DC <: PSY.DCSource,
    P <: PSY.FrequencyEstimator,
    F <: PSY.Filter,
}
    return
end

function _update_inner_vars!(
    device_states::AbstractArray{<:ACCEPTED_REAL_TYPES},
    ::AbstractArray{<:ACCEPTED_REAL_TYPES},
    ω_sys::ACCEPTED_REAL_TYPES,
    inner_vars::AbstractArray{<:ACCEPTED_REAL_TYPES},
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
    θ = atan(V_I, V_R)

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
    Z_source_sq = R_source^2 + X_source^2

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
            return Vr_cnv, Vi_cnv
        else
            return Vr_inv, Vi_inv
        end
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
        #Compute aux currents
        I_aux_r = (R_source * Vr_cnv + X_source * Vi_cnv) / Z_source_sq
        I_aux_i = (R_source * Vi_cnv - X_source * Vr_cnv) / Z_source_sq
        Ir_filt = Ir_cnv - I_aux_r
        Ii_filt = Ii_cnv - I_aux_i
    end

    #Update inner_vars
    inner_vars[Ir_cnv_var] = Ir_cnv
    inner_vars[Ii_cnv_var] = Ii_cnv
    inner_vars[Vr_cnv_var] = Vr_cnv
    inner_vars[Vi_cnv_var] = Vi_cnv
    inner_vars[Vr_filter_var] = Vr_cnv
    inner_vars[Vi_filter_var] = Vi_cnv
    inner_vars[Ir_filter_var] = Ir_cnv
    inner_vars[Ii_filter_var] = Ii_cnv
    inner_vars[Ir_inv_var] = Ir_filt
    inner_vars[Ii_inv_var] = Ii_filt
    return
end

function _update_inner_vars!(
    device_states::AbstractArray{<:ACCEPTED_REAL_TYPES},
    ::AbstractArray{<:ACCEPTED_REAL_TYPES},
    ω_sys::ACCEPTED_REAL_TYPES,
    inner_vars::AbstractArray{<:ACCEPTED_REAL_TYPES},
    dynamic_device::DynamicWrapper{
        PSY.DynamicInverter{
            PSY.RenewableEnergyVoltageConverterTypeA,
            O,
            IC,
            DC,
            P,
            PSY.LCLFilter,
        },
    },
) where {
    O <: PSY.OuterControl,
    IC <: PSY.InnerControl,
    DC <: PSY.DCSource,
    P <: PSY.FrequencyEstimator,
}
    filter_ix = get_local_state_ix(dynamic_device, PSY.LCLFilter)
    filter_states = @view device_states[filter_ix]
    Ir_cnv = filter_states[1]
    Ii_cnv = filter_states[2]
    Vr_filter = filter_states[3]
    Vi_filter = filter_states[4]
    Ir_filter = filter_states[5]
    Ii_filter = filter_states[6]

    #Update inner_vars
    inner_vars[Ir_cnv_var] = Ir_cnv
    inner_vars[Ii_cnv_var] = Ii_cnv
    inner_vars[Vr_filter_var] = Vr_filter
    inner_vars[Vi_filter_var] = Vi_filter
    inner_vars[Ir_inv_var] = Ir_filter
    inner_vars[Ii_inv_var] = Ii_filter
    inner_vars[Ir_filter_var] = Ir_filter
    inner_vars[Ii_filter_var] = Ii_filter
    return
end

function device_mass_matrix_entries!(
    mass_matrix::AbstractArray,
    dynamic_device::DynamicWrapper{T},
) where {
    T <: Union{PSY.SingleCageInductionMachine, PSY.SimplifiedSingleCageInductionMachine},
}
    global_index = get_global_index(dynamic_device)
    mass_matrix_induction_entries!(mass_matrix, dynamic_device, global_index)
    return
end

function mass_matrix_induction_entries!(
    mass_matrix,
    ind::DynamicWrapper{T},
    global_index::ImmutableDict{Symbol, Int64},
) where {
    T <: Union{PSY.SingleCageInductionMachine, PSY.SimplifiedSingleCageInductionMachine},
}
    @debug "Using default mass matrix entries $ind"
end

"""
Model of 5-state (SingleCageInductionMachine) induction motor in Julia.
Refer to "Analysis of Electric Machinery and Drive Systems" by Paul Krause,
Oleg Wasynczuk and Scott Sudhoff for the equations
"""
function device!(
    device_states::AbstractArray{T},
    output_ode::AbstractArray{T},
    voltage_r::T,
    voltage_i::T,
    current_r::AbstractArray{T},
    current_i::AbstractArray{T},
    global_vars::AbstractArray{T},
    ::AbstractArray{T},
    dynamic_wrapper::DynamicWrapper{PSY.SingleCageInductionMachine},
    t,
) where {T <: ACCEPTED_REAL_TYPES}
    Sbase = get_system_base_power(dynamic_wrapper)
    f0 = get_system_base_frequency(dynamic_wrapper)
    if get_connection_status(dynamic_wrapper) < 1.0
        output_ode .= zero(T)
        return
    end

    # get speed of system's reference frame
    sys_ω = global_vars[GLOBAL_VAR_SYS_FREQ_INDEX]

    # get states
    ψ_qs = device_states[1]
    ψ_ds = device_states[2]
    ψ_qr = device_states[3]
    ψ_dr = device_states[4]
    ωr = device_states[5]

    #Get parameters
    dynamic_device = get_device(dynamic_wrapper)
    R_s = PSY.get_R_s(dynamic_device)
    X_ls = PSY.get_X_ls(dynamic_device)
    R_r = PSY.get_R_r(dynamic_device)
    X_lr = PSY.get_X_lr(dynamic_device)
    A = PSY.get_A(dynamic_device)
    B = PSY.get_B(dynamic_device)
    C = PSY.get_C(dynamic_device)
    H = PSY.get_H(dynamic_device)
    base_power = PSY.get_base_power(dynamic_device)
    B_sh = PSY.get_B_shunt(dynamic_device)
    τ_m0 = PSY.get_τ_ref(dynamic_device)
    X_ad = PSY.get_X_ad(dynamic_device)
    X_aq = PSY.get_X_aq(dynamic_device)

    # voltages in QD
    v_qs = voltage_i
    v_ds = voltage_r
    v_qr = zero(T)
    v_dr = zero(T)

    #Additional Fluxes
    ψ_mq = X_aq * (ψ_qs / X_ls + ψ_qr / X_lr) # (4.14-15) in Krause
    ψ_md = X_ad * (ψ_ds / X_ls + ψ_dr / X_lr) # (4.14-16) in Krause

    # Stator motor currents in QD
    i_qs = 1 / X_ls * (ψ_qs - ψ_mq) # (4.14-1) in Krause
    i_ds = 1 / X_ls * (ψ_ds - ψ_md) # (4.14-2) in Krause

    # Electric Torque
    τ_e = ψ_ds * i_qs - ψ_qs * i_ds    # (4.14-18) in Krause

    #Compute ODEs
    output_ode[1] = 2.0 * pi * f0 * (v_qs - sys_ω * ψ_ds - R_s * i_qs)  # (4.14-9) in Krause
    output_ode[2] = 2.0 * pi * f0 * (v_ds + sys_ω * ψ_qs - R_s * i_ds) # (4.14-10) in Krause
    output_ode[3] =
        2.0 * pi * f0 * (v_qr - (sys_ω - ωr) * ψ_dr + R_r / X_lr * (ψ_mq - ψ_qr)) # (4.14-12) in Krause
    output_ode[4] =
        2.0 * pi * f0 * (v_dr + (sys_ω - ωr) * ψ_qr + R_r / X_lr * (ψ_md - ψ_dr)) # (4.14-13) in Krause
    output_ode[5] = 1.0 / (2.0 * H) * (τ_e - τ_m0 * (A * ωr^2 + B * ωr + C)) # (4.14-19) in Krause and (7.39) in Kundur

    #Update current
    current_r[1] -= (base_power / Sbase) * (i_ds - v_qs * B_sh)  # in system base
    current_i[1] -= (base_power / Sbase) * (i_qs + v_ds * B_sh)  # in system base
    return
end

"""
Model of 3-state (SimplifiedSingleCageInductionMachine) induction motor in Julia.
Based on the 3rd order model derived in Prabha Kundur's Book and the
equations in "Analysis of Electric Machinery and Drive Systems" by Paul Krause,
Oleg Wasynczuk and Scott Sudhoff.
"""
function device!(
    device_states::AbstractArray{T},
    output_ode::AbstractArray{T},
    voltage_r::T,
    voltage_i::T,
    current_r::AbstractArray{T},
    current_i::AbstractArray{T},
    global_vars::AbstractArray{T},
    ::AbstractArray{T},
    dynamic_wrapper::DynamicWrapper{PSY.SimplifiedSingleCageInductionMachine},
    t,
) where {T <: ACCEPTED_REAL_TYPES}
    Sbase = get_system_base_power(dynamic_wrapper)
    f0 = get_system_base_frequency(dynamic_wrapper)
    if get_connection_status(dynamic_wrapper) < 1.0
        output_ode .= zero(T)
        return
    end

    # get speed of system's reference frame
    sys_ω = global_vars[GLOBAL_VAR_SYS_FREQ_INDEX]

    # get states
    ψ_qr = device_states[1]
    ψ_dr = device_states[2]
    ωr = device_states[3]

    #Get parameters
    dynamic_device = get_device(dynamic_wrapper)
    R_s = PSY.get_R_s(dynamic_device)
    X_m = PSY.get_X_m(dynamic_device)
    R_r = PSY.get_R_r(dynamic_device)
    A = PSY.get_A(dynamic_device)
    B = PSY.get_B(dynamic_device)
    C = PSY.get_C(dynamic_device)
    H = PSY.get_H(dynamic_device)
    base_power = PSY.get_base_power(dynamic_device)
    B_sh = PSY.get_B_shunt(dynamic_device)
    τ_m0 = PSY.get_τ_ref(dynamic_device)
    X_rr = PSY.get_X_rr(dynamic_device)
    X_p = PSY.get_X_p(dynamic_device)

    # voltages in QD
    v_qs = voltage_i
    v_ds = voltage_r
    v_qr = zero(T)
    v_dr = zero(T)

    # Stator and rotor currents in QD
    i_qs =
        1 / (R_s^2 + (sys_ω * X_p)^2) * (
            (R_s * v_qs - sys_ω * X_p * v_ds) -
            (R_s * sys_ω * X_m / X_rr * ψ_dr + sys_ω * X_p * sys_ω * X_m / X_rr * ψ_qr)
        )
    i_ds =
        1 / (R_s^2 + (sys_ω * X_p)^2) * (
            (R_s * v_ds + sys_ω * X_p * v_qs) -
            (-R_s * sys_ω * X_m / X_rr * ψ_qr + sys_ω * X_p * sys_ω * X_m / X_rr * ψ_dr)
        )
    i_qr = (ψ_qr - X_m * i_qs) / X_rr # derived from 4th row of (4.5.37) in Krause
    i_dr = (ψ_dr - X_m * i_ds) / X_rr # # derived from 5th row of (4.5.37) in Krause

    # Electric Torque
    τ_e = ψ_qr * i_dr - ψ_dr * i_qr     # (4.8-6) in Krause

    #Compute ODEs
    output_ode[1] = 2.0 * pi * f0 * (v_qr - (sys_ω - ωr) * ψ_dr - R_r * i_qr) # (4.5-25) in Krause
    output_ode[2] = 2.0 * pi * f0 * (v_dr + (sys_ω - ωr) * ψ_qr - R_r * i_dr) # (4.5-26) in Krause
    output_ode[3] = 1.0 / (2.0 * H) * (τ_e - τ_m0 * (A * ωr^2 + B * ωr + C)) # (4.14-19) in Krause and (7.39) in Kundur

    #Update current
    current_r[1] -= (base_power / Sbase) * (i_ds - v_qs * B_sh)  # in system base
    current_i[1] -= (base_power / Sbase) * (i_qs + v_ds * B_sh)  # in system base
    return
end

function device_mass_matrix_entries!(
    mass_matrix::AbstractArray,
    dynamic_device::DynamicWrapper{PSY.CSVGN1},
)
    global_index = get_global_index(dynamic_device)
    mass_matrix_csvgn1_entries!(mass_matrix, dynamic_device, global_index)
    return
end

function mass_matrix_csvgn1_entries!(
    mass_matrix,
    csvgn1::DynamicWrapper{PSY.CSVGN1},
    global_index::ImmutableDict{Symbol, Int64},
)
    @debug "Using default mass matrix entries $csvgn1"
end

"""
Model of Static Shunt Compensator: CSVGN1.
"""
function device!(
    device_states::AbstractArray{T},
    output_ode::AbstractArray{T},
    voltage_r::T,
    voltage_i::T,
    current_r::AbstractArray{T},
    current_i::AbstractArray{T},
    global_vars::AbstractArray{T},
    ::AbstractArray{T},
    dynamic_wrapper::DynamicWrapper{PSY.CSVGN1},
    t,
) where {T <: ACCEPTED_REAL_TYPES}
    Sbase = get_system_base_power(dynamic_wrapper)
    V_ref = get_V_ref(dynamic_wrapper)
    # TODO: V_abs is the voltage magnitude on the high side of generator step-up transformer, if present.
    V_abs = sqrt(voltage_r^2 + voltage_i^2)

    if get_connection_status(dynamic_wrapper) < 1.0
        output_ode .= zero(T)
        return
    end

    # get states
    thy = device_states[1]
    vr1 = device_states[2]
    vr2 = device_states[3]

    #Get parameters
    dynamic_device = get_device(dynamic_wrapper)
    K = PSY.get_K(dynamic_device)
    T1 = PSY.get_T1(dynamic_device)
    T2 = PSY.get_T2(dynamic_device)
    T3 = PSY.get_T3(dynamic_device)
    T4 = PSY.get_T4(dynamic_device)
    T5 = PSY.get_T5(dynamic_device)
    Rmin = PSY.get_Rmin(dynamic_device)
    Vmax = PSY.get_Vmax(dynamic_device)
    Vmin = PSY.get_Vmin(dynamic_device)
    Cbase = PSY.get_CBase(dynamic_device)
    # FIXME: base_power is changed to system's base_power when a CSVGN1 is attached to a Source using add_component!()
    # Temporarily, to avoid that, set_dynamic_injector!() could be used
    Mbase = PSY.get_base_power(dynamic_device)
    Rbase = Mbase

    # Regulator
    T3_ll = T1 + T2 # time parameter for the lead-lag block
    T4_ll = T1 * T2 # time parameter for the lead-lag block
    T1_ll = T3 + T4 # time parameter for the lead-lag block
    T2_ll = T3 * T4 # time parameter for the lead-lag block

    v_ll, dvr1_dt, dvr2_dt = lead_lag_2nd_nonwindup(
        K * (V_abs - V_ref),
        vr1,
        vr2,
        T1_ll,
        T2_ll,
        T3_ll,
        T4_ll,
        Vmin,
        Vmax,
        )

    # Thyristor
    y_r, dthy_dt = low_pass_nonwindup(v_ll, thy, 1.0, T5, Rmin / Rbase, 1.0)

    # Admittance output
    Y = Cbase / Sbase - y_r * Mbase / Sbase

    #Compute ODEs
    output_ode[1] = dthy_dt
    output_ode[2] = dvr1_dt
    output_ode[3] = dvr2_dt

    #Update current
    current_r[1] = (Mbase / Sbase) * Y * voltage_i # in system base
    current_i[1] = -(Mbase / Sbase) * Y * voltage_r # in system base
    return
end

function device_mass_matrix_entries!(
    mass_matrix::AbstractArray,
    dynamic_device::DynamicWrapper{PSY.ActiveConstantPowerLoad},
)
    global_index = get_global_index(dynamic_device)
    device = get_device(dynamic_device)
    bool_mm_value = PSY.get_is_filter_differential(device)
    f0 = get_system_base_frequency(dynamic_device)
    ωb = 2 * pi * f0
    mass_matrix[global_index[:ir_cnv], global_index[:ir_cnv]] =
        bool_mm_value * PSY.get_lf(device) / ωb
    mass_matrix[global_index[:ii_cnv], global_index[:ii_cnv]] =
        bool_mm_value * PSY.get_lf(device) / ωb
    mass_matrix[global_index[:vr_filter], global_index[:vr_filter]] =
        bool_mm_value * PSY.get_cf(device) / ωb
    mass_matrix[global_index[:vi_filter], global_index[:vi_filter]] =
        bool_mm_value * PSY.get_cf(device) / ωb
    mass_matrix[global_index[:ir_filter], global_index[:ir_filter]] =
        bool_mm_value * PSY.get_lg(device) / ωb
    mass_matrix[global_index[:ii_filter], global_index[:ii_filter]] =
        bool_mm_value * PSY.get_lg(device) / ωb
    return
end

"""
Model of 12-state Active Constant Power Load in Julia.
Based on the paper `Malicious Control of an Active Load in an Islanded Mixed-Source Microgrid`
by C. Roberts, U. Markovic, D. Arnold and D. Callaway.
"""
function device!(
    device_states::AbstractArray{T},
    output_ode::AbstractArray{T},
    voltage_r::T,
    voltage_i::T,
    current_r::AbstractArray{T},
    current_i::AbstractArray{T},
    global_vars::AbstractArray{T},
    ::AbstractArray{T},
    dynamic_wrapper::DynamicWrapper{PSY.ActiveConstantPowerLoad},
    t,
) where {T <: ACCEPTED_REAL_TYPES}
    Sbase = get_system_base_power(dynamic_wrapper)
    f0 = get_system_base_frequency(dynamic_wrapper)
    V_ref = get_V_ref(dynamic_wrapper)
    if get_connection_status(dynamic_wrapper) < 1.0
        output_ode .= zero(T)
        return
    end

    # get speed of system's reference frame
    sys_ω = global_vars[GLOBAL_VAR_SYS_FREQ_INDEX]

    # get states
    θ_pll = device_states[1]
    ϵ_pll = device_states[2]
    η = device_states[3]
    v_dc = device_states[4]
    γd = device_states[5]
    γq = device_states[6]
    Ir_cnv = device_states[7]
    Ii_cnv = device_states[8]
    Vr_filter = device_states[9]
    Vi_filter = device_states[10]
    Ir_filter = device_states[11]
    Ii_filter = device_states[12]

    # Voltages in Network Ref Frame
    V_tR = voltage_r
    V_tI = voltage_i

    # Transformation to Ref Frame
    I_dq_cnv = ri_dq(θ_pll + pi / 2) * [Ir_cnv; Ii_cnv]

    #Get parameters
    dynamic_device = get_device(dynamic_wrapper)
    r_load = PSY.get_r_load(dynamic_device)
    c_dc = PSY.get_c_dc(dynamic_device)
    rf = PSY.get_rf(dynamic_device)
    lf = PSY.get_lf(dynamic_device)
    cf = PSY.get_cf(dynamic_device)
    rg = PSY.get_rg(dynamic_device)
    lg = PSY.get_lg(dynamic_device)
    kp_pll = PSY.get_kp_pll(dynamic_device)
    ki_pll = PSY.get_ki_pll(dynamic_device)
    kpv = PSY.get_kpv(dynamic_device)
    kiv = PSY.get_kiv(dynamic_device)
    kpc = PSY.get_kpc(dynamic_device)
    kic = PSY.get_kic(dynamic_device)
    base_power = PSY.get_base_power(dynamic_device)

    # Compute PLL expressions
    V_dq_pll = ri_dq(θ_pll + pi / 2) * [Vr_filter; Vi_filter]
    Δω_pi, dϵ_dt = pi_block(V_dq_pll[2], ϵ_pll, kp_pll, ki_pll)
    ω_pll = Δω_pi + sys_ω

    # Compute DC side output
    Id_ref, dη_dt = pi_block(V_ref - v_dc, η, kpv, kiv)
    Iq_ref = get_Q_ref(dynamic_wrapper)

    # Compute AC controller expressions
    Vd_ref_uncomp, dγd_dt = pi_block(-Id_ref + I_dq_cnv[d], γd, kpc, kic)
    Vq_ref_uncomp, dγq_dt = pi_block(-Iq_ref + I_dq_cnv[q], γq, kpc, kic)
    Vd_cnv = Vd_ref_uncomp + ω_pll * lf * I_dq_cnv[q]
    Vq_cnv = Vq_ref_uncomp - ω_pll * lf * I_dq_cnv[d]

    Vr_cnv, Vi_cnv = dq_ri(θ_pll + pi / 2) * [Vd_cnv; Vq_cnv]

    # Converter Power
    P_cnv = Vr_cnv * Ir_cnv + Vi_cnv * Ii_cnv

    #Compute ODEs

    ## PLL Dynamics
    output_ode[1] = 2 * pi * f0 * Δω_pi
    output_ode[2] = dϵ_dt
    ## DC Voltage Controller
    output_ode[3] = dη_dt
    output_ode[4] = (2 * pi * f0 / c_dc) * (P_cnv / v_dc - v_dc / r_load)
    ## AC Current Controller
    output_ode[5] = dγd_dt
    output_ode[6] = dγq_dt
    ## AC Dynamics
    #𝜕id_c/𝜕t
    output_ode[7] = (Vr_filter - Vr_cnv - rf * Ir_cnv + lf * sys_ω * Ii_cnv)
    #𝜕iq_c/𝜕t
    output_ode[8] = (Vi_filter - Vi_cnv - rf * Ii_cnv - lf * sys_ω * Ir_cnv)
    #LCL Capacitor (internal state)
    #𝜕vd_o/𝜕t
    output_ode[9] = (Ir_filter - Ir_cnv + cf * sys_ω * Vi_filter)
    #𝜕vq_o/𝜕t
    output_ode[10] = (Ii_filter - Ii_cnv - cf * sys_ω * Vr_filter)
    #Grid Inductance (internal state)
    #𝜕id_o/𝜕t
    output_ode[11] = (V_tR - Vr_filter - rg * Ir_filter + lg * sys_ω * Ii_filter)
    #𝜕iq_o/𝜕t
    output_ode[12] = (V_tI - Vi_filter - rg * Ii_filter - lg * sys_ω * Ir_filter)

    #Update current
    current_r[1] -= Ir_filter * base_power / Sbase  # in system base
    current_i[1] -= Ii_filter * base_power / Sbase  # in system base
    return
end

function device_mass_matrix_entries!(
    mass_matrix::AbstractArray,
    dynamic_device::DynamicWrapper{PSY.AggregateDistributedGenerationA},
)
    global_index = get_global_index(dynamic_device)
    mass_matrix_dera_entries!(mass_matrix, dynamic_device, global_index)
    return
end

function mass_matrix_dera_entries!(
    mass_matrix,
    dera::DynamicWrapper{PSY.AggregateDistributedGenerationA},
    global_index::ImmutableDict{Symbol, Int64},
)
    ddera = get_device(dera)
    Freq_Flag = PSY.get_Freq_Flag(get_device(dera))
    if Freq_Flag == 1
        mass_matrix[global_index[:Vmeas], global_index[:Vmeas]] = PSY.get_T_rv(ddera)
        mass_matrix[global_index[:Pmeas], global_index[:Pmeas]] = PSY.get_Tp(ddera)
        mass_matrix[global_index[:Q_V], global_index[:Q_V]] = PSY.get_T_iq(ddera)
        mass_matrix[global_index[:Mult], global_index[:Mult]] = PSY.get_Tv(ddera)
        mass_matrix[global_index[:Fmeas], global_index[:Fmeas]] = PSY.get_Trf(ddera)
        mass_matrix[global_index[:Pord], global_index[:Pord]] = PSY.get_Tpord(ddera)
    else
        mass_matrix[global_index[:Vmeas], global_index[:Vmeas]] = PSY.get_T_rv(ddera)
        mass_matrix[global_index[:Pmeas], global_index[:Pmeas]] = PSY.get_Tp(ddera)
        mass_matrix[global_index[:Q_V], global_index[:Q_V]] = PSY.get_T_iq(ddera)
        mass_matrix[global_index[:Mult], global_index[:Mult]] = PSY.get_Tv(ddera)
        mass_matrix[global_index[:Fmeas], global_index[:Fmeas]] = PSY.get_Trf(ddera)
    end
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
    dynamic_device::DynamicWrapper{PSY.AggregateDistributedGenerationA},
    t,
) where {T <: ACCEPTED_REAL_TYPES}
    Freq_Flag = PSY.get_Freq_Flag(get_device(dynamic_device))
    _mdl_ode_AggregateDistributedGenerationA!(
        device_states,
        output_ode,
        Val(Freq_Flag),
        voltage_r,
        voltage_i,
        current_r,
        current_i,
        global_vars,
        inner_vars,
        dynamic_device,
        t,
    )
    return
end

#####################################################
### Auxiliary ODE calculations via Flags dispatch ###
#####################################################

#Freq_Flag = 0
function _mdl_ode_AggregateDistributedGenerationA!(
    device_states::AbstractArray{T},
    output_ode::AbstractArray{T},
    ::Val{0},
    voltage_r::T,
    voltage_i::T,
    current_r::AbstractArray{T},
    current_i::AbstractArray{T},
    global_vars::AbstractArray{T},
    inner_vars::AbstractArray{T},
    dynamic_device::DynamicWrapper{PSY.AggregateDistributedGenerationA},
    t,
) where {T <: ACCEPTED_REAL_TYPES}
    sys_ω = global_vars[GLOBAL_VAR_SYS_FREQ_INDEX]
    Sbase = get_system_base_power(dynamic_device)
    Vt = sqrt(voltage_r^2 + voltage_i^2)

    #Obtain References (from wrapper and device)
    Pfa_ref = PSY.get_Pfa_ref(get_device(dynamic_device))
    P_ref = get_P_ref(dynamic_device)
    Q_ref = get_Q_ref(dynamic_device)
    V_ref = get_V_ref(dynamic_device)

    #Get flags
    Pf_Flag = PSY.get_Pf_Flag(get_device(dynamic_device))

    #Get device states
    Vmeas = device_states[1]
    Pmeas = device_states[2]
    Q_V = device_states[3]
    Iq = device_states[4]
    Mult = device_states[5]
    Fmeas = device_states[6]
    Ip = device_states[7]

    Ip_cmd = Ip
    Iq_cmd = Iq

    #Get parameters
    T_rv = PSY.get_T_rv(get_device(dynamic_device))
    Trf = PSY.get_Trf(get_device(dynamic_device))
    (dbd1, dbd2) = PSY.get_dbd_pnts(get_device(dynamic_device))
    K_qv = PSY.get_K_qv(get_device(dynamic_device))
    Tp = PSY.get_Tp(get_device(dynamic_device))
    T_iq = PSY.get_T_iq(get_device(dynamic_device))

    Tg = PSY.get_Tg(get_device(dynamic_device))
    rrpwr = PSY.get_rrpwr(get_device(dynamic_device))
    Tv = PSY.get_Tv(get_device(dynamic_device))
    Iq_lim = PSY.get_Iq_lim(get_device(dynamic_device))
    basepower = PSY.get_base_power(get_device(dynamic_device))
    base_power_ratio = basepower / Sbase

    #STATE Vmeas
    _, dVmeas_dt = low_pass_mass_matrix(Vt, Vmeas, 1.0, T_rv)
    #STATE Q_V
    if Pf_Flag == 1
        _, dQ_V_dt =
            low_pass_mass_matrix(tan(Pfa_ref) * Pmeas / max(Vmeas, 0.01), Q_V, 1.0, T_iq)
    elseif Pf_Flag == 0
        _, dQ_V_dt = low_pass_mass_matrix(Q_ref / max(Vmeas, 0.01), Q_V, 1.0, T_iq)
    else
        @error @error "Unsupported value of PQ_Flag"
    end

    #STATE Iq
    Ip_min, Ip_max, Iq_min, Iq_max =
        current_limit_logic(get_device(dynamic_device), Ip_cmd, Iq_cmd)
    Iq_input =
        clamp(
            clamp(
                deadband_function(V_ref - Vmeas, dbd1, dbd2) * K_qv,
                Iq_lim[:min],
                Iq_lim[:max],
            ) + Q_V,
            Iq_min,
            Iq_max,
        ) * Mult
    _, dIq_dt = low_pass(Iq_input, Iq, 1.0, Tg)

    VMult = 1.0
    FMult = 1.0
    _, dMult_dt = low_pass_mass_matrix(VMult * FMult, Mult, 1.0, Tv)

    #STATE Fmeas
    _, dFmeas_dt = low_pass_mass_matrix(sys_ω, Fmeas, 1.0, Trf)

    if Ip >= 0
        Rup = abs(rrpwr)
        Rdown = -Inf
    else
        Rdown = -abs(rrpwr)
        Rup = Inf
    end

    #STATE Ip
    Ip_input = clamp(P_ref / max(Vmeas, 0.01), Ip_min, Ip_max) * Mult
    Ip_limited, dIp_dt =
        low_pass_nonwindup_ramp_limits(Ip_input, Ip, 1.0, Tg, -Inf, Inf, Rdown, Rup)

    #STATE Pmeas
    _, dPmeas_dt = low_pass_mass_matrix(P_ref, Pmeas, 1.0, Tp)

    #Update ODEs
    output_ode[1] = dVmeas_dt
    output_ode[2] = dPmeas_dt
    output_ode[3] = dQ_V_dt
    output_ode[4] = dIq_dt
    output_ode[5] = dMult_dt
    output_ode[6] = dFmeas_dt
    output_ode[7] = dIp_dt

    #Calculate output current
    θ = atan(voltage_i / voltage_r)
    Iq_neg = -Iq
    I_r = real(complex(Ip_limited, Iq_neg) * exp(im * θ))
    I_i = imag(complex(Ip_limited, Iq_neg) * exp(im * θ))
    current_r[1] = I_r * base_power_ratio
    current_i[1] = I_i * base_power_ratio
end

#Freq_Flag = 1
function _mdl_ode_AggregateDistributedGenerationA!(
    device_states::AbstractArray{T},
    output_ode::AbstractArray{T},
    ::Val{1},
    voltage_r::T,
    voltage_i::T,
    current_r::AbstractArray{T},
    current_i::AbstractArray{T},
    global_vars::AbstractArray{T},
    inner_vars::AbstractArray{T},
    dynamic_device::DynamicWrapper{PSY.AggregateDistributedGenerationA},
    t,
) where {T <: ACCEPTED_REAL_TYPES}
    sys_ω = global_vars[GLOBAL_VAR_SYS_FREQ_INDEX]
    Sbase = get_system_base_power(dynamic_device)
    Vt = sqrt(voltage_r^2 + voltage_i^2)

    #Obtain References (from wrapper and device)
    Pfa_ref = PSY.get_Pfa_ref(get_device(dynamic_device))
    P_ref = get_P_ref(dynamic_device)
    Q_ref = get_Q_ref(dynamic_device)
    V_ref = get_V_ref(dynamic_device)
    ω_ref = get_ω_ref(dynamic_device)

    #Get flags
    Pf_Flag = PSY.get_Pf_Flag(get_device(dynamic_device))

    #Get device states
    Vmeas = device_states[1]
    Pmeas = device_states[2]
    Q_V = device_states[3]
    Iq = device_states[4]
    Mult = device_states[5]
    Fmeas = device_states[6]
    PowerPI = device_states[7]
    dPord = device_states[8]
    Pord = device_states[9]
    Ip = device_states[10]

    Ip_cmd = Ip
    Iq_cmd = Iq

    #Get parameters
    T_rv = PSY.get_T_rv(get_device(dynamic_device))
    Trf = PSY.get_Trf(get_device(dynamic_device))
    (dbd1, dbd2) = PSY.get_dbd_pnts(get_device(dynamic_device))
    K_qv = PSY.get_K_qv(get_device(dynamic_device))
    Tp = PSY.get_Tp(get_device(dynamic_device))
    T_iq = PSY.get_T_iq(get_device(dynamic_device))
    D_dn = PSY.get_D_dn(get_device(dynamic_device))
    D_up = PSY.get_D_up(get_device(dynamic_device))
    (fdbd1, fdbd2) = PSY.get_fdbd_pnts(get_device(dynamic_device))
    fe_lim = PSY.get_fe_lim(get_device(dynamic_device))
    P_lim = PSY.get_P_lim(get_device(dynamic_device))
    dP_lim = PSY.get_dP_lim(get_device(dynamic_device))
    Tpord = PSY.get_Tpord(get_device(dynamic_device))
    Kpg = PSY.get_Kpg(get_device(dynamic_device))
    Kig = PSY.get_Kig(get_device(dynamic_device))

    Tg = PSY.get_Tg(get_device(dynamic_device))
    rrpwr = PSY.get_rrpwr(get_device(dynamic_device))
    Tv = PSY.get_Tv(get_device(dynamic_device))
    Iq_lim = PSY.get_Iq_lim(get_device(dynamic_device))
    basepower = PSY.get_base_power(get_device(dynamic_device))
    base_power_ratio = basepower / Sbase

    #STATE Vmeas
    _, dVmeas_dt = low_pass_mass_matrix(Vt, Vmeas, 1.0, T_rv)
    #STATE Q_V
    if Pf_Flag == 1
        _, dQ_V_dt =
            low_pass_mass_matrix(tan(Pfa_ref) * Pmeas / max(Vmeas, 0.01), Q_V, 1.0, T_iq)
    elseif Pf_Flag == 0
        _, dQ_V_dt = low_pass_mass_matrix(Q_ref / max(Vmeas, 0.01), Q_V, 1.0, T_iq)
    else
        @error "Unsupported value of PQ_Flag"
    end

    #STATE Iq
    Ip_min, Ip_max, Iq_min, Iq_max =
        current_limit_logic(get_device(dynamic_device), Ip_cmd, Iq_cmd)
    Iq_input =
        clamp(
            clamp(
                deadband_function(V_ref - Vmeas, dbd1, dbd2) * K_qv,
                Iq_lim[:min],
                Iq_lim[:max],
            ) + Q_V,
            Iq_min,
            Iq_max,
        ) * Mult
    _, dIq_dt = low_pass(Iq_input, Iq, 1.0, Tg)

    VMult = 1.0
    FMult = 1.0
    _, dMult_dt = low_pass_mass_matrix(VMult * FMult, Mult, 1.0, Tv)

    #STATE Fmeas
    _, dFmeas_dt = low_pass_mass_matrix(sys_ω, Fmeas, 1.0, Trf)

    if Ip >= 0
        Rup = abs(rrpwr)
        Rdown = -Inf
    else
        Rdown = -abs(rrpwr)
        Rup = Inf
    end

    #STATE Pmeas
    _, dPmeas_dt = low_pass_mass_matrix(P_ref, Pmeas, 1.0, Tp)

    #STATE PowerPI
    PowerPI_input = clamp(
        min(deadband_function(ω_ref - Fmeas, fdbd1, fdbd2) * D_dn, 0.0) +
        max(deadband_function(ω_ref - Fmeas, fdbd1, fdbd2) * D_up, 0.0) - Pmeas + P_ref,
        fe_lim[:min],
        fe_lim[:max],
    )
    _, dPowerPI_dt =
        pi_block_nonwindup(PowerPI_input, PowerPI, Kpg, Kig, P_lim[:min], P_lim[:max])

    #STATE dPord
    if dPowerPI_dt > dP_lim[:max]
        ddPord_dt = dP_lim[:max]
    elseif dPowerPI_dt < dP_lim[:min]
        ddPord_dt = dP_lim[:min]
    else
        ddPord_dt = dPowerPI_dt
    end

    #State Pord
    Pord_limited, dPord_dt =
        low_pass_nonwindup_mass_matrix(dPord, Pord, 1.0, Tpord, P_lim[:min], P_lim[:max])

    #STATE Ip
    Ip_input = clamp(Pord_limited / max(Vmeas, 0.01), Ip_min, Ip_max) * Mult
    Ip_limited, dIp_dt =
        low_pass_nonwindup_ramp_limits(Ip_input, Ip, 1.0, Tg, -Inf, Inf, Rdown, Rup)

    #Update ODEs
    output_ode[1] = dVmeas_dt
    output_ode[2] = dPmeas_dt
    output_ode[3] = dQ_V_dt
    output_ode[4] = dIq_dt
    output_ode[5] = dMult_dt
    output_ode[6] = dFmeas_dt
    output_ode[7] = dPowerPI_dt
    output_ode[8] = ddPord_dt
    output_ode[9] = dPord_dt
    output_ode[10] = dIp_dt

    #Calculate output current
    θ = atan(voltage_i / voltage_r)
    Iq_neg = -Iq
    I_r = real(complex(Ip_limited, Iq_neg) * exp(im * θ))
    I_i = imag(complex(Ip_limited, Iq_neg) * exp(im * θ))
    current_r[1] = I_r * base_power_ratio
    current_i[1] = I_i * base_power_ratio
end
