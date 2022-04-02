function mass_matrix_converter_entries!(
    mass_matrix,
    converter::C,
    global_index::Base.ImmutableDict{Symbol, Int64},
) where {C <: PSY.Converter}
    @debug "Using default mass matrix entries $C"
end

function mdl_converter_ode!(
    device_states::AbstractArray{<:ACCEPTED_REAL_TYPES},
    output_ode::AbstractArray{<:ACCEPTED_REAL_TYPES},
    inner_vars::AbstractArray{<:ACCEPTED_REAL_TYPES},
    dynamic_device::DynamicWrapper{
        PSY.DynamicInverter{PSY.AverageConverter, O, IC, DC, P, F},
    },
) where {
    O <: PSY.OuterControl,
    IC <: PSY.InnerControl,
    DC <: PSY.DCSource,
    P <: PSY.FrequencyEstimator,
    F <: PSY.Filter,
}

    #Obtain inner variables for component
    md = inner_vars[md_var]
    mq = inner_vars[mq_var]
    Vdc = inner_vars[Vdc_var]
    θ_oc = inner_vars[θ_oc_var]

    #Transform reference frame to grid reference frame
    m_ri = dq_ri(θ_oc + pi / 2) * [md; mq]

    #Update inner_vars
    inner_vars[Vr_cnv_var] = m_ri[R] * Vdc
    inner_vars[Vi_cnv_var] = m_ri[I] * Vdc
    return
end

function mdl_converter_ode!(
    device_states::AbstractArray{<:ACCEPTED_REAL_TYPES},
    output_ode::AbstractArray{<:ACCEPTED_REAL_TYPES},
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
    #Obtain inner variables for component
    V_R = inner_vars[Vr_inv_var]
    V_I = inner_vars[Vi_inv_var]
    V_t = sqrt(V_R^2 + V_I^2)
    θ = atan(V_I, V_R)
    Ip_cmd = inner_vars[Id_ic_var]
    Iq_cmd = inner_vars[Iq_ic_var]

    #Get Converter parameters
    converter = PSY.get_converter(dynamic_device)
    T_g = PSY.get_T_g(converter)
    Rrpwr = PSY.get_Rrpwr(converter)
    Brkpt = PSY.get_Brkpt(converter)
    Zerox = PSY.get_Zerox(converter)
    Lvpl1 = PSY.get_Lvpl1(converter)
    Vo_lim = PSY.get_Vo_lim(converter)
    Lv_pnt0, Lv_pnt1 = PSY.get_Lv_pnts(converter)
    # Io_lim = PSY.get_Io_lim(converter)
    T_fltr = PSY.get_T_fltr(converter)
    K_hv = PSY.get_K_hv(converter)
    Iqr_min, Iqr_max = PSY.get_Iqr_lims(converter)
    Lvpl_sw = PSY.get_Lvpl_sw(converter)
    Q_ref = PSY.get_Q_ref(converter)
    R_source = PSY.get_R_source(converter)
    X_source = PSY.get_X_source(converter)

    #Obtain indices for component w/r to device
    local_ix = get_local_state_ix(dynamic_device, PSY.RenewableEnergyConverterTypeA)
    #Define internal states for Converter
    internal_states = @view device_states[local_ix]
    Ip = internal_states[1]
    Iq = internal_states[2]
    Vmeas = internal_states[3]

    # Compute additional variables
    # Active Power Part
    #Update Ip ramp limits
    Rp_up = Ip >= 0 ? Rrpwr : Inf
    Rp_dn = Ip >= 0 ? -Inf : -Rrpwr
    #Saturate Ip if LVPL is active
    Ip_sat = Ip
    Ip_binary = 1.0
    if Lvpl_sw == 1
        LVPL = get_LVPL_gain(Vmeas, Zerox, Brkpt, Lvpl1)
        Ip_sat = Ip <= LVPL ? Ip : LVPL
        Ip_binary = Ip <= LVPL ? 1.0 : 0.0
    end
    Ip_in = clamp(Ip_cmd - Ip_sat, Rp_dn, Rp_up)
    #Get Low Voltage Active Current Management Gain
    G_lv = get_LV_current_gain(V_t, Lv_pnt0, Lv_pnt1)

    # Reactive Power Part
    #Update Iq ramp limits
    Rq_up = Q_ref >= 0 ? Iqr_max : Inf
    Rq_dn = Q_ref >= 0 ? -Inf : Iqr_min
    Iq_in = clamp(Iq_cmd - Iq, Rq_dn, Rq_up)
    Iq_extra = max(K_hv * (V_t - Vo_lim), 0.0)
    Id_cnv = G_lv * Ip_sat
    Iq_cnv = -Iq - Iq_extra
    #Reference Transformation
    Ir_cnv = Id_cnv * cos(θ) - Iq_cnv * sin(θ)
    Ii_cnv = Id_cnv * sin(θ) + Iq_cnv * cos(θ)

    ###
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
            lf_ratio = lf / Zf_mag_sq
            Vr_cnv =
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

    #Compute converter voltage
    Vr_cnv, Vi_cnv = V_cnv_calc(Ir_cnv, Ii_cnv, V_R, V_I)

    #Update ODEs
    output_ode[local_ix[1]] = Ip_binary * (1.0 / T_g) * Ip_in #(Ip_cmd - Ip)
    output_ode[local_ix[2]] = (1.0 / T_g) * Iq_in # (Iq_cmd - Iq)
    output_ode[local_ix[3]] = (1.0 / T_fltr) * (V_t - Vmeas)

    #Update inner_vars
    inner_vars[Ir_cnv_var] = Ir_cnv
    inner_vars[Ii_cnv_var] = Ii_cnv
    inner_vars[Vr_cnv_var] = Vr_cnv
    inner_vars[Vi_cnv_var] = Vi_cnv
    return
end
